using Random
using LightGraphs

################ This files contains tools for the setups of the formulations, heuristics, decompositions, and separation algorithm



function lecture_graphe(filename) #to obtain an adjacency matrix
    
    open(filename) do f
        h = split(read(filename,String)," ")
        g = fill(0,1,1)
        if(h[1] == "c" || h[1] =="p")
        
           for (i,line) in enumerate(eachline(f))
              x = split(line," ")
              if(x[1]=="p")
                  n = parse(Int,x[3])
                  g = fill(0,n,n)
              elseif(x[1] == "e")
                  v_1 = parse(Int, x[2])
                  v_2 = parse(Int, x[3])
                  g[v_1,v_2] = 1
                  g[v_2,v_1] = 1
               end
           end
        end
    
        return g
    end
    
    
    
end


function cleanGraph2plex(g,infbound)  ## Remove as a preprocessing the vertices that cannot belong to the optimal 2-plex as their neighborhood is too small


	bool = false
	V = Array{Int}(undef,0)
	for i = 1 : nv(g)
	
		if(size( neighbors(g,i))[1] >= infbound-2)
		 push!(V,i)
		else
		
			bool = true
		end
	
	end
	
	if(bool)
		return cleanGraph2plex(induced_subgraph(g,V)[1],infbound)
	else
		return g
	
	end

	
end

function cleanGraphco2plex(g,infbound)
	
	
	#return complement(cleanGraph2plex(complement(g),infbound))
	
	bool = false
	V = Array{Int}(undef,0)
	for i = 1 : nv(g)
	
		if(nv(g) - size( neighbors(g,i))[1] > infbound-2)
		 push!(V,i)
		else
		
			bool = true
		end
	
	end
	
	if(bool)
		return cleanGraphco2plex(induced_subgraph(g,V)[1],infbound)
	else
		return g
	
	end


	
end

function iscandidate(g,i,K)
	if(size(K)[1]==0)
		return true
	end
	m= 0
	if(i in K )
		return false
	
	end

	for j in K
		if(!has_edge(g,i,j))
			m = m+1
			
			
			
			
			if(m == 2)
				return false
			end
			
			if(size([k for k in neighbors(g,j) if k in K])[1] < size(K)[1]-1 )
				return false
			end
		end

	end
	return true
end

function getBiggestLowestN(g1,list) ## returns the maximal and minimal size of a neighborhood
	min = nv(g1)
	max = 0

	
	for v =1 :  size(list)[1]
	
		s = size(neighbors(g1,v))[1]
		if(s < min )
			min = s
	
		end
		if(s>max)
			max = s
		end
	end
	return max,min
end

function heuristic_max_2_plex(g,alpha,K)# builds a maximum 2-plex of a graph according 
	
	candidates = [i for i = 1 : nv(g)]
	
	for i = 1 : nv(g)
	
		candidates = [j for j in candidates if (iscandidate(g, j , K))]
		
		if(size(candidates)[1] == 0)
		
			break
		end
		g1,list = induced_subgraph(g,candidates)
		max, min = getBiggestLowestN(g1,list)
		threshold = min + alpha*(max - min)
		candidates_serious = [list[j] for j =1 :  size(list)[1] if size(neighbors(g1,j))[1] >= threshold]
		push!(K, rand(candidates_serious))
		
		
	
	end
	


	K = localSearch2Plex(g,K)  ##
	# is_2plex(g,K)
	return K
end

function is_2plex(g,k) ## tests wether a set of vertices is a 2-plex
	println("test")
	g1 = induced_subgraph(g,k)[1]

	for i = 1 : nv(g1)
	
		if(size(neighbors(g1,i))[1] < size(k)[1]-2)
			println("not 2 plex")
			return false
		end
	
	
	
	end
	return true

end
function localSearch2Plex(g,K) 
	
	notInK = [i for i = 1 : nv(g) if (! (i in K))]
	
	
	
	for i in K
		K1 = [j for j in K if j != i]

		for j in notInK
			
			if(iscandidate(g,j,K1))
				K2 = K1
				push!(K2,j)
				for k in notInK
					if(iscandidate(g,k,K2))
						push!(K2,k)
						return localSearch2Plex(g,K2)
					
					end		
				
				end
			end
			
		end
	
	end
	return K
end

function heuristic_max_clique(g,costs)

	return heuristic_max_stable(complement(g),costs)

end

function heuristic_max_stable(g,costs)  # computes a maximum sized stable set


	couples = [(i,costs[i]) for i = 1 : nv(g)]
	sort!(couples, by = x-> x[2],rev= true)
	
	result = []
	cost = 0
	while(size(couples)[1] != 0)
		push!(result, couples[1][1])
		cost = cost + couples[1][2]
		for i in neighbors(g,couples[1][1])
			
			for j = 2 : size(couples)[1]
				if (i==couples[j][1])
					deleteat!(couples,j)
					break
				
				end
						
			end

		end
		deleteat!(couples,1)		
	
	end	

	
	return result,cost

end  

function completeIntoMaxStable(g, stable)

	vertices = []
	
	for i = 1 : nv(g)
		if(! (i in stable))
		
			independent = true
			for j in stable
				if(has_edge(g, i, j)) 
					independent = false
				end
		
			end
			if(independent)
				push!(stable, i)
			
			end
			
		end
	
	
	end
	
		
		
		

	
	
	return stable
end

function getNeighborhoods(g)
	result= []
	neighborhoods = []

	for i =1 : nv(g)
		r= Array{Int}(undef, 0)
		for j in neighbors(g,i)
			if(!(j in r))
				push!(r,j)
			end
			for k in neighbors(g,j)
				if(k!= i && (!( k in r)))
					push!(r, k)
				end
			end
		
		end
		push!(r,i)
		push!(result,induced_subgraph(g,r))
	
	
	end
	
	return result
	

end



function coloringVertices(tg_bar,g,index,bb)
	#edges = Random.shuffle([i for i = 1 : ne(g)])
	result=  []
	
	for e in edges(g)
		stable=[src(e),dst(e),index[src(e),dst(e)]]
		if(bb)
		for i in vertices(g)
			
			if(has_edge(g,src(e),i)  && i!=dst(e) && i != src(e))
				push!(stable,index[src(e),i])
			end
			if(has_edge(g,dst(e),i)  && i!=src(e) && i!= dst(e))
				push!(stable,index[dst(e),i])
			end
		end
		end

		#println("edge " * string(u) * " " * string(v))
		#stable = completeIntoMaxStable(tg_bar,stable)
		push!(result,stable)

	
	end
		
	return result
end

function getStableOfSubset(g,set)
	result = []
	for i in set	
		independent = true
		for j in result
			if(has_edge(g,j,i))
				independent = false
			end
		end
		if(independent)
			push!(result, i)
		end	
	end
	return result
end


function augmentedTotalGraph(g) ### Computes the utter graph

	result = SimpleGraph(nv(g)+ne(g))
	
	index = zeros(Int64, nv(g), nv(g))
	compteur = 1+nv(g)
	for uv in edges(g)
		if(src(uv) != dst(uv) )
		index[src(uv),dst(uv)] = compteur
		index[dst(uv),src(uv)] = compteur
		
		add_edge!(result, dst(uv),src(uv)) 
		add_edge!(result,src(uv),compteur) 
		add_edge!(result,dst(uv),compteur)
		compteur = compteur+1
		end
	end
	
	compteur = compteur-1
	#total graph
	
	for v in vertices(g)
		for u in neighbors(g,v)
			if(u!= v)
				for i in neighbors(g,u)
					if(i != v && u !=i && !has_edge(result,v,index[u,i]))

							add_edge!(result,v,index[u,i])
				
						
				
					end
				end
			end
		end
		
	
	end
	# il faut ajouter les arêtes entre sommets arrêtes
	
	for e in edges(g)
		for i = 1 : nv(g)
			if(i != src(e) && i != dst(e))		
		
				if(has_edge(g,src(e),i) && ! has_edge(result,index[src(e),i], index[src(e),dst(e)]))
					add_edge!(result,index[src(e),i], index[src(e),dst(e)] )
				elseif( has_edge(g,dst(e),i)  && !has_edge(result,index[dst(e),i], index[src(e),dst(e)]))
					add_edge!(result,index[dst(e),i], index[src(e),dst(e)] )
				
				end
			end
		
		end
	
	
	end
	
	for e1 in edges(g)
		for e2 in edges(g)
		
			if(e1 != e2)
				if(has_edge(g,src(e1),src(e2)) && !has_edge(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)]))
					add_edge!(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)])
				elseif(has_edge(g,src(e1),dst(e2)) && !has_edge(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)]))
					add_edge!(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)])			
				elseif(has_edge(g,dst(e1),dst(e2)) && !has_edge(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)]))
					add_edge!(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)])
				elseif(has_edge(g,dst(e1),src(e2))&& !has_edge(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)]))
					add_edge!(result,index[src(e1),dst(e1)],index[src(e2),dst(e2)])
				
				end
			end
		end	
	end	
	
	
	
	
	return result, index

end
