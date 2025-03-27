using CPLEX
using JuMP, Cbc
using MathOptInterface
using Dates

include("graph.jl")


####################### This file contains the differents models described in the corresponding papers
####################### Each function correspond to a different model, some are reinforced with primal constraints: the separation is declared inside the function
####################### Each function set a timeout of one hour to which the time spent solving the previous bricks of the decomposition is substracted
####################### Each function returns the statistics  of the resolution

function SimpleModel(g1,costfunction,separate,totalTime) ############ compact model of Balasundaram

	g = g1
	bar_g = complement(g)

	
	
	result = Model( CPLEX.Optimizer)
	MOI.set(result, MOI.Silent(), true)
	set_time_limit_sec(result, max(1,3600.0-totalTime))  ## Set a timeout depending on how much have been spent on the other bricks of the instance's 'decomposition
	@variable(result,   1 >= x[1 : nv(g)]>= 0,Int)  # One variable by vertex

	for i = 1 : nv(g) # Each constraint tells that if a vertex is taken in the solution then at most one of its neighbors can be taken too
		@constraint(result, sum(x[j] for j in neighbors(g,i) ) + (size(neighbors(g,i))[1] -1 ) *x[i] <= size(neighbors(g,i))[1])
	end
	nb_separate = 0
	function separatingCliqueIq(cb_data) ## Dynamic separation of  clique inequalities, 
	        callback_called = true
	        x_vals = callback_value.(Ref(cb_data), x)
        	
        	
        	clique, value = heuristic_max_stable(bar_g,x_vals)
        	
		if value > 2
			nb_separate = nb_separate+1
		    con = @build_constraint(
		        sum(x[i] for i in clique ) <= 2
		    )
		    MOI.submit(result, MOI.UserCut(cb_data), con)
		end
	end
	
	
	if(separate)  ### If separation should be added
		MOI.set(result, MOI.UserCutCallback(), separatingCliqueIq)
	end
	
	
	#MOI.set(result, MOI.HeuristicCallback(), heuristicc)
	@objective(result ,Max, sum(x[i]*costfunction[i] for i = 1 : nv(g)))  ## Set objective coefficient
	
	
	
	optimize!(result)


	return termination_status(result), solve_time(result), node_count(result),nb_separate,objective_value(result) ## return statistics 

end





function StableSetInTotalNoSep(g1,costfunction,with_edges,totalTime) ## Computes a maximum stable set of the utter graph
	
	g = g1
	bar_g = complement(g)
	
	tg,index = augmentedTotalGraph(g)
	tg_bar = complement(tg)		

	result = Model( CPLEX.Optimizer)
	MOI.set(result, MOI.Silent(), true)

	@variable(result,   1 >= z[1 : nv(tg)]>= 0,Int)
	tt = now()
	coloringTG = coloringVertices(tg_bar,g,index,true)  # a clique cover of TG so the formulation is valid
	tt = (now()-tt).value / 1000
	
	set_time_limit_sec(result, 3600.0-totalTime-tt) ## Set time limit 
	
	for i in coloringTG  ## Add sufficiently enough clique constraints
		@constraint(result, sum(z[j] for j in i) <= 1)
	end
	
	if(with_edges) ## Should we add constraints for each edge of the utter graph
		for e in edges(tg)
			@constraint(result, z[src(e)] + z[dst(e)] <= 1)
		end
	end
	
	## Vertices of the utter graph associated with edges of the input graph cost double
	@objective(result ,Max, sum(z[i]*costfunction[i] for i = 1 : nv(g)) + sum(z[index[src(e),dst(e)]]*(costfunction[src(e)] + costfunction[dst(e)])  for e in edges(g)))	
		
	optimize!(result)
	return termination_status(result), solve_time(result)+ tt, node_count(result),0,objective_value(result)  ## return statistics

end

function StableSetInTotal(g1,costfunction,with_edges,totalTime) ## Computes a stable set of the utter graph reinforced by clique inequalities

	
	g = g1
	bar_g = complement(g)
	nb_separate = 0
	tg,index = augmentedTotalGraph(g)
	tg_bar = complement(tg)	

	result = Model( CPLEX.Optimizer)
	MOI.set(result, MOI.Silent(), true)
		@variable(result,   1 >= z[1 : nv(tg)]>= 0,Int)
	tt = 0	
	


	if(with_edges)
		for e in edges(tg)
			@constraint(result, z[src(e)] + z[dst(e)] <= 1)
		end
		set_time_limit_sec(result, max(1,3600.0-totalTime))
	else

		tt = now()
		coloringTG = coloringVertices(tg_bar,g,index,true)  # a clique cover of TG
		tt= (now()-tt).value / 1000	
		set_time_limit_sec(result, max(1,3600.0-totalTime-tt))
		for i in coloringTG
			@constraint(result, sum(z[j] for j in i) <= 1)
		end

	end
	
	
	@objective(result ,Max, sum(z[i]*costfunction[i] for i = 1 : nv(g)) + sum(z[index[src(e),dst(e)]]*(costfunction[src(e)] + costfunction[dst(e)])   for e in edges(g) ))
	
	function separatingCliqueIq(cb_data) ## The clique separation algorithm
	        callback_called = true
	        z_vals = callback_value.(Ref(cb_data), z)
        	
        	clique, value = heuristic_max_stable(tg_bar,z_vals)
        	
		if value > 1
	 	    nb_separate = nb_separate +1
		    con = @build_constraint(
		        sum(z[i] for i in clique) <= 1
		    )
		    MOI.submit(result, MOI.UserCut(cb_data), con)
		end
	end
	if(with_edges)
		MOI.set(result, MOI.UserCutCallback(), separatingCliqueIq)
	end

	optimize!(result)
	return termination_status(result), solve_time(result)+tt, node_count(result),nb_separate,objective_value(result)

end

function coeffcientOfEdge(g,index, nb,stable) #returns the coefficient of a variable in the clique constraint.
	e = collect(edges(g))[nb-nv(g)]
	if(src(e) in stable)
		if(dst(e) in stable)
			return -1
		else
			return 0
		end
	else
		if(dst(e) in stable)
			return 0
		else
			return 1
		end
	end
end

function ProjectedForm(g1,costfunction,separate,totalTime) ## The projected formulation  with a separation of clique inequalities

	
	g = g1
	bar_g = complement(g)
	tg,index = augmentedTotalGraph(g)
	tg_bar = complement(tg)	
	
	index_edge = []
	for e in edges(g)
		push!(index_edge, [src(e),dst(e)])
	
	end
	
	result = Model( CPLEX.Optimizer)
	MOI.set(result, MOI.Silent(), true)
	
	@variable(result,   1 >= x[1 : nv(g)]>= 0,Int)
	@variable(result,   1 >= y[1 : ne(g)]>= 0,Int)
	tt = 0
	if(separate)
		coloringTG = coloringVertices(tg_bar,g,index,false)
		set_time_limit_sec(result, max(1,3600.0-totalTime))
	else
		tt = now()
		coloringTG = coloringVertices(tg_bar,g,index,true)
		tt = (now() - tt).value / 1000
		set_time_limit_sec(result, max(1,3600.0-totalTime-tt))


	end

	
	for i in coloringTG #clique constraints
		@constraint(result, sum(x[j] for j in i if j <= nv(g)) + sum(y[j-nv(g)]*coeffcientOfEdge(g,index,j,i) for j in i if j >nv(g) ) <= 1 )
	end
	
	for i = 1 : nv(g) #neighborhood constraints
		@constraint(result, sum(y[index[j,i]-nv(g)] for j in neighbors(g,i) ) <= x[i])	
	end

	
	
	@objective(result ,Max, sum(x[j] * costfunction[j] for j = 1 : nv(g)) ) 
	nb_separate = 0
	deg = degree(g)
	

	function separatingCliqueIq(cb_data)
	        callback_called = true
	        x_vals = callback_value.(Ref(cb_data), x)
        	y_vals = callback_value.(Ref(cb_data), y)
       		
        	z_vals = [x_vals[i] -  debile(deg, i ,g,y_vals,index) for i = 1:nv(g)]
        	co = vcat(z_vals,y_vals)
        	clique, value = heuristic_max_stable(tg_bar,co)
        	
		if value > 1
			nb_separate = nb_separate+1
		    con = @build_constraint(
		        sum(x[i] for i in clique if i<=nv(g)) - sum(y[i-nv(g)] *coeffcientOfEdge(g,index,i,clique)  for i in clique if (i>nv(g)))  <= 1
		    )
		    MOI.submit(result, MOI.UserCut(cb_data), con)
		end
	end
	
	if(separate)
		MOI.set(result, MOI.UserCutCallback(), separatingCliqueIq)
	end
	
	optimize!(result)
	return termination_status(result), solve_time(result)+tt, node_count(result),nb_separate,objective_value(result)
   
end

function debile(deg, i ,g,y_vals,index)
	if(deg[i]>0)
		return sum(y_vals[index[i,j]-nv(g)] for j in neighbors(g,i))

	else
		return 0
	
	end


end

function SimpleModelSepNeigh(g1,costfunction,separate,totalTime) ## Implementation of the simple model reinforced with nbeighborhood inequalities separation
	
g = g1
	bar_g = complement(g)

	
	
	result = Model( CPLEX.Optimizer)
	MOI.set(result, MOI.Silent(), true)
	set_time_limit_sec(result, max(1,3600.0-totalTime))
	set_attribute(result, "CPXPARAM_Simplex_DynamicRows", 1)
	@variable(result,   1 >= x[1 : nv(g)]>= 0,Int)

	for i = 1 : nv(g)
		@constraint(result, sum(x[j] for j in neighbors(g,i) ) + (size(neighbors(g,i))[1] -1 ) *x[i] <= size(neighbors(g,i))[1])
	end
	nb_separate = 0
	
	
	function separatingNeighIq(cb_data)

	        callback_called = true
	        x_vals = callback_value.(Ref(cb_data), x)
	        
        	
        	for u = 1: nv(g) # for each vertex compute the inequality that is the most violated
        		if(size(neighbors(g,u))[1] > 1)
		    		sorted_neigh = sort([ [v, x_vals[v]]  for v in neighbors(g,u) if v!=u], by = x -> x[2])
		    		set = [sorted_neigh[end][1], sorted_neigh[end-1][1]]
		    		sum = x_vals[u] + sorted_neigh[end][2] + sorted_neigh[end-1][2]
		    		
		    		for k = length(sorted_neigh)-2:-1:1
		    			if sorted_neigh[k][2] + x_vals[u] > 1
		    				sum = sum + sorted_neigh[k][2] + x_vals[u]
		    				push!(set,sorted_neigh[k][1])
		    			
		    			else 
		    				break
		    			end
		    		
		    		end
		    		
		    		if(sum > length(set)+ 1e-6)

						nb_separate = nb_separate+1
						con = @build_constraint(sum(x[Int(i)] for i in set ) + (length(set)-1) * x[u] <= length(set))
						
						MOI.submit(result, MOI.UserCut(cb_data), con)
						
					
					end
        		end
        	
        		
        	end
        	return 

	end
	
	
	MOI.set(result, MOI.UserCutCallback(), separatingNeighIq)

	
	@objective(result ,Max, sum(x[i]*costfunction[i] for i = 1 : nv(g)))
	
	
	
	optimize!(result)


	return termination_status(result), solve_time(result), node_count(result),nb_separate,objective_value(result)
end

	
	



	

