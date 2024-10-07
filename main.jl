using CPLEX
using JuMP, Cbc
using LightGraphs
using Random
using Shuffle
using MathOptInterface


include("graph.jl")
include("model_ilp.jl")


Random.seed!(3213)




function solveCo_2_Plex(fun::Function,g,costs,bool,name,pseudo_partition)
	

	bar_g = complement(g)

	
	alg = String(Symbol(fun))*"_"*string(bool)
	
	totalSolved=0
	totalTime = 0
	totalNode = 0
	totalSeparate = 0
	max_opt = 0
	println("start solving: " * String(Symbol(fun)) * " " * name)
	for i in pseudo_partition
		if(size(i[2])[1] != 0)
			solved, time, nb_node,nb_separate, opt = fun(complement(i[1]),[costs[j] for j in i[2]],bool,totalTime)
			if(Int(solved) == 1) #TerminationStatusCode.OPTIMAL)
				totalSolved=totalSolved+1
			end
			totalTime = totalTime+time
			totalNode = totalNode+nb_node
			totalSeparate = totalSeparate + nb_separate
			if(max_opt < opt)
			
				max_opt = opt
			end
			if(totalTime > 3600)
				file = open(alg*".txt","a")
				write(file, name *"&" * string(totalSolved == nv(g))*"&"*string(totalTime)*"&"*string(totalNode)*"&"*string(totalSeparate)*"&"* string(max_opt) * "\n")
				close(file)
				return
			end
		end
	end
		

		
	file = open(alg*".txt","a")
	write(file, name *"&" * string(totalSolved == nv(g))*"&"*string(totalTime)*"&"*string(totalNode)*"&"*string(totalSeparate)*"&"* string(max_opt) * "\n")
	close(file)
	
end

function solveInstanceSimple(fName, co_2_plex)
	
	open(fName) do f
		for (i,line) in enumerate(eachline(f))
		
			if(co_2_plex)
				g = SimpleGraph(lecture_graphe("data/"*line))
				
			else
				g = complement(SimpleGraph(lecture_graphe("data/"*line)))
			
			end
			
			bar_g = complement(g)
	
			ss=size(heuristic_max_2_plex(bar_g, 0.7,[]))[1]
			g = cleanGraphco2plex(g,ss)
			g1 = copy(g)
			bar_g1 = complement(g1)
                	pseudo_partition = getNeighborhoods(bar_g1)
			pseudo_partition1 = []
			for gg in pseudo_partition
				bar_gg = complement(gg[1])
				ss=size(heuristic_max_2_plex(bar_gg, 0.7,[]))[1]
				subg = cleanGraphco2plex(gg[1],ss)
				push!(pseudo_partition1, [subg, gg[2]])
			end

			costs = [1 for i = 1 : nv(g)]
	       		g1=copy(g)
			solveCo_2_Plex(SimpleModel,g1,costs,true,line,pseudo_partition1)
			
			
		end
        end
end

function solveInstanceSimpleNeigh(fName, co_2_plex)
	
	open(fName) do f
		for (i,line) in enumerate(eachline(f))
		
			if(co_2_plex)
				g = SimpleGraph(lecture_graphe("data/"*line))
				
			else
				g = complement(SimpleGraph(lecture_graphe("data/"*line)))
			
			end
			
			bar_g = complement(g)
	
			ss=size(heuristic_max_2_plex(bar_g, 0.7,[]))[1]
			g = cleanGraphco2plex(g,ss)
			g1 = copy(g)
			bar_g1 = complement(g1)
                	pseudo_partition = getNeighborhoods(bar_g1)
			pseudo_partition1 = []
			for gg in pseudo_partition
				bar_gg = complement(gg[1])
				ss=size(heuristic_max_2_plex(bar_gg, 0.7,[]))[1]
				subg = cleanGraphco2plex(gg[1],ss)
				push!(pseudo_partition1, [subg, gg[2]])
			end

			costs = [1 for i = 1 : nv(g)]
	       		g1=copy(g)
			solveCo_2_Plex(SimpleModelSepNeigh,g1,costs,true,line,pseudo_partition1)
			
			
		end
        end
end

function solveInstanceStableTrue(fName, co_2_plex)
	
	open(fName) do f
		for (i,line) in enumerate(eachline(f))
		
			if(co_2_plex)
				g = SimpleGraph(lecture_graphe("data/"*line))
				
			else
				g = complement(SimpleGraph(lecture_graphe("data/"*line)))
			
			end
			
			bar_g = complement(g)
	
			ss=size(heuristic_max_2_plex(bar_g, 0.7,[]))[1]
			g = cleanGraphco2plex(g,ss)
			g1 = copy(g)
			bar_g1 = complement(g1)
                	pseudo_partition = getNeighborhoods(bar_g1)
			pseudo_partition1 = []
			for gg in pseudo_partition
				bar_gg = complement(gg[1])
				ss=size(heuristic_max_2_plex(bar_gg, 0.7,[]))[1]
				subg = cleanGraphco2plex(gg[1],ss)
				push!(pseudo_partition1, [subg, gg[2]])
			end

			costs = [1 for i = 1 : nv(g)]
			g1 = copy(g)
			solveCo_2_Plex(StableSetInTotal,g1,costs,true,line,pseudo_partition1)
			
		end
        end
end
function solveInstanceStableFalse(fName, co_2_plex)
	
	open(fName) do f
		for (i,line) in enumerate(eachline(f))
		
			if(co_2_plex)
				g = SimpleGraph(lecture_graphe("data/"*line))
				
			else
				g = complement(SimpleGraph(lecture_graphe("data/"*line)))
			
			end
			
			bar_g = complement(g)
	
			ss=size(heuristic_max_2_plex(bar_g, 0.7,[]))[1]
			g = cleanGraphco2plex(g,ss)
			g1 = copy(g)
			bar_g1 = complement(g1)
                	pseudo_partition = getNeighborhoods(bar_g1)
			pseudo_partition1 = []
			for gg in pseudo_partition
				bar_gg = complement(gg[1])
				ss=size(heuristic_max_2_plex(bar_gg, 0.7,[]))[1]
				subg = cleanGraphco2plex(gg[1],ss)
				push!(pseudo_partition1, [subg, gg[2]])
			end

			costs = [1 for i = 1 : nv(g)]
	       		g1 = copy(g)
			solveCo_2_Plex(StableSetInTotal,g1,costs,false,line,pseudo_partition1)
			
		end
        end
end
function solveInstanceProjectedFormTrue(fName, co_2_plex)
	
	open(fName) do f
		for (i,line) in enumerate(eachline(f))
		
			if(co_2_plex)
				g = SimpleGraph(lecture_graphe("data/"*line))
				
			else
				g = complement(SimpleGraph(lecture_graphe("data/"*line)))
			
			end
			
			bar_g = complement(g)
	
			ss=size(heuristic_max_2_plex(bar_g, 0.7,[]))[1]
			g = cleanGraphco2plex(g,ss)
			g1 = copy(g)
			bar_g1 = complement(g1)
                	pseudo_partition = getNeighborhoods(bar_g1)
			pseudo_partition1 = []
			for gg in pseudo_partition
				bar_gg = complement(gg[1])
				ss=size(heuristic_max_2_plex(bar_gg, 0.7,[]))[1]
				subg = cleanGraphco2plex(gg[1],ss)
				push!(pseudo_partition1, [subg, gg[2]])
			end

			costs = [1 for i = 1 : nv(g)]
			g1 = copy(g)
			
			solveCo_2_Plex(ProjectedForm,g1,costs,true,line,pseudo_partition1)
			

		end
        end
end
function solveInstanceProjectedFormFalse(fName, co_2_plex)
	
	open(fName) do f
		for (i,line) in enumerate(eachline(f))
		
			if(co_2_plex)
				g = SimpleGraph(lecture_graphe("data/"*line))
				
			else
				g = complement(SimpleGraph(lecture_graphe("data/"*line)))
			
			end
			
			bar_g = complement(g)
	
			ss=size(heuristic_max_2_plex(bar_g, 0.7,[]))[1]
			g = cleanGraphco2plex(g,ss)
			g1 = copy(g)
			bar_g1 = complement(g1)
                	pseudo_partition = getNeighborhoods(bar_g1)
			pseudo_partition1 = []
			for gg in pseudo_partition
				bar_gg = complement(gg[1])
				ss=size(heuristic_max_2_plex(bar_gg, 0.7,[]))[1]
				subg = cleanGraphco2plex(gg[1],ss)
				push!(pseudo_partition1, [subg, gg[2]])
			end

			costs = [1 for i = 1 : nv(g)]
	       		g1 = copy(g)
			solveCo_2_Plex(ProjectedForm,g1,costs,false,line,pseudo_partition1)

		end
        end
end



#writeChiDOIRel(df_doi,df_basic)
#y=getNeighborhoods(g)


#SimpleModel(g,costs,true)
#SimpleModel(g,costs,false)

g = SimpleGraph(lecture_graphe("data/miles250.col"))
costs = [1 for i = 1 :nv(g)]

ProjectedForm(g,costs,true,0.0)




#g = SimpleGraph(lecture_graphe("data/1-FullIns_5.col"))
#costs = [1 for i = 1 :nv(g)]

#println("c lu")
#solveCo_2_Plex(StableSetInTotalNoSep,g,costs, true,"yo")
#g1 = SimpleGraph(lecture_graphe("data/1-FullIns_3.col"))

#println(SimpleModelSepNeigh(g1,[1 for i in vertices(g1)],false,0))
#println(SimpleModel(g1,[1 for i in vertices(g1)],false,0))
#solveInstanceSimple("instances.txt",true)
solveInstanceSimpleNeigh("instances.txt",true)


#solveInstanceStableTrue("instancesStableTrue.txt",true)

#solveInstanceStableFalse("instancesStableFalse.txt",true)

#solveInstanceProjectedFormTrue("instancesProjectedFormTrue.txt",true)

#solveInstanceProjectedFormFalse("instancesProjectedFormFalse.txt",true)
