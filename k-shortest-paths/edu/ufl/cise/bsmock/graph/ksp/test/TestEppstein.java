package edu.ufl.cise.bsmock.graph.ksp.test;
import java.lang.Math;


import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import edu.ufl.cise.bsmock.graph.Graph;
import edu.ufl.cise.bsmock.graph.ksp.Eppstein;
import edu.ufl.cise.bsmock.graph.util.Path;

import java.util.HashMap;
import java.util.Map;
import java.util.Hashtable;
import java.util.List;

/**
 * Test of Eppstein's algorithm for computing the K shortest paths between two nodes in a graph.
 *
 * Copyright (C) 2015  Brandon Smock (dr.brandon.smock@gmail.com, GitHub: bsmock)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created by Brandon Smock on October 7, 2015.
 * Last updated by Brandon Smock on December 24, 2015.
 * 

 */
public class TestEppstein {

    public static void main(String args[]) throws IOException {
        /* Uncomment any of these example tests */
    	String run_mode = "DEPLOY";

    	if (run_mode.equals("TEST")) {
			String graphFilename_w, graphFilename_uw, outFilename_uw,outFilename_w, sourceNodes, targetNode,source_contr;
			int K;
			double tolerance;
			int len_threshold;
			
			// *********************************** Testing code ****************************************
			String g_type = "unweighted";
			
			graphFilename_uw = "edu/ufl/cise/bsmock/graph/ksp/test/all_same_weight_input-graph-ss-genemaniaplus-a0_01.txt";
			outFilename_uw = "edu/ufl/cise/bsmock/graph/ksp/test/out_uw.txt";
			
			graphFilename_w = "edu/ufl/cise/bsmock/graph/ksp/test/shortest-path-input-graph-ss-genemaniaplus-a0_01.txt";
			outFilename_w = "edu/ufl/cise/bsmock/graph/ksp/test/out_w.txt";
			
			sourceNodes = "[9161]";
			targetNode = "6162";
			K = 1000;
			len_threshold = 4;
			
			if (g_type == "unweighted") {
				shortest_paths_by_length_on_unweighted(graphFilename_uw,outFilename_uw,sourceNodes,targetNode,K, len_threshold);
			}
			else if  (g_type == "weighted") {
				shortest_paths_by_length_on_weighted(graphFilename_w, outFilename_w, sourceNodes,  targetNode, K ,len_threshold, outFilename_uw);
			}
    	}
      
    	else if (run_mode.equals("DEPLOY")) {
	    	// *********************************** Integratable Portion of the code ******************************
	    	String threshold_type = args[0];
			if (threshold_type.equals("pathnum"))
			{
				int i = 1;
				String graphFilename = args[i];
				String outFilename = args[i + 1];
				String sourceNodes = args[i + 2];  //comma separated list of nodes e.g. "["source_1", "source_2"]"
				String targetNode = args[i + 3];
				Integer K = Integer.parseInt(args[i + 4]);
				shortest_paths_by_pathnum(graphFilename, outFilename,
						sourceNodes, targetNode, K);
			}
			else if (threshold_type.equals("pathlen"))
			{
				String graph_type = args[1];
				//	        System.out.println(graph_type);
				if (graph_type.equals("unweighted")) {
					int i = 2;
					String graphFilename = args[i];
					String outFilename = args[i + 1];
					String sourceNodes = args[i + 2];  //comma separated list of nodes e.g. "["source_1", "source_2"]"
					String targetNode = args[i + 3];
					Integer K = Integer.parseInt(args[i + 4]);
					Integer len_threshold = Integer.parseInt(args[i + 5]);
					shortest_paths_by_length_on_unweighted(graphFilename, outFilename, sourceNodes, targetNode, K, len_threshold);
					//		        System.out.println("Finished unweighted");

				} else if (graph_type.equals("weighted")) {
					int i = 2;
					String graphFilename = args[i];
					String outFilename = args[i + 1];
					String sourceNodes = args[i + 2];
					String targetNode = args[i + 3];
					Integer K = Integer.parseInt(args[i + 4]);
					Integer len_threshold = Integer.parseInt(args[i + 5]);
					String len_threshold_filename = args[i + 6];

					shortest_paths_by_length_on_weighted(graphFilename, outFilename,
							sourceNodes, targetNode, K, len_threshold, len_threshold_filename);
					//		        System.out.println("Finished weighted");

				}
			}
    	}
//  
    }

	public static void shortest_paths_by_pathnum(String graphFilename, String outFilename,
																String source, String target, int k) throws IOException
	{

		/* Read graph from file */
		Graph graph = new Graph(graphFilename);

		Eppstein eppsteinAlgorithm = new Eppstein();

//      source is a string with multiple comma separated sources
		source = source.replace("[", "");
		source = source.replace("]", "");
		String[] sources = source.split("," ,-2);
//        System.out.println("Sources: "+ sources[0]);
		PrintWriter outWriterAll = new PrintWriter(new FileOutputStream(new File(outFilename), false));
		List<Path> ksp;

		int source_count = 0;
		while (source_count<sources.length) {
			HashMap<Integer, Integer> paths_per_len = new HashMap<Integer, Integer>();
			String s = sources[source_count];
			s = s.replace(" ","");
//        	System.out.println("s: "+ s);
			ksp = eppsteinAlgorithm.ksp(graph, s, target, k);

			for (Path p : ksp) {
				System.out.println(p);
				outWriterAll.append(s +" "+target +" "+ p.toString()+"\n");
			}
			source_count+=1;
		}

		ksp=null;
		graph=null;
		outWriterAll.close();

	}

	public static void shortest_paths_by_tolerance(String graphFilename, String outFilename,
    		String source, String target, int k, String source_contr, double tolerance) throws IOException
    {
        /* Read graph from file */
        Graph graph = new Graph(graphFilename);

        Eppstein eppsteinAlgorithm = new Eppstein();
        
//      source is a string with multiple comma separated sources
        source = source.replace("[", "");
        source = source.replace("]", "");
        String[] sources = source.split("," ,-2);
        
        source_contr = source_contr.replace("[", "");
        source_contr = source_contr.replace("]", "");
        String[] source_contrs = source_contr.split("," ,-2);
//        System.out.print(sources);
        
        /* Output the K shortest paths */
//        PrintWriter outWriter = new PrintWriter(tempOutFilename);
        PrintWriter outWriterAll = new PrintWriter(new FileOutputStream(new File(outFilename), false));
    	List<Path> ksp;
    	
    	int source_count = 0;
    	int cur_k = k;
    	double sum=0;
    	double prev_sum=0;
        while (source_count<sources.length) {
        	String s = sources[source_count];
        	s = s.replace(" ","");
//        	System.out.println("s: "+ s);
        	ksp = eppsteinAlgorithm.ksp(graph, s, target, cur_k);
        	
        	//check if the fraction of contributions from shortest paths sum up to (1-tolerance). If not then run ksp with higher K.
        	for (Path p : ksp) {
        		sum+=Math.pow(10, -p.getTotalCost());
        	}
//        	
//        	System.out.println(sum);
//        	System.out.println(prev_sum);
//        	System.out.println(sum-prev_sum);

//        	System.out.println(source_contr);
//        	System.out.println(sum/source_contr);
        	
        	double contr = Double.parseDouble(source_contrs[source_count]);
        	if ((sum/contr)>=(1-tolerance)| ((sum-prev_sum)<Math.pow(10,-10))){ // at some point if adding new path does not change the contribution much, then also stop.
        		
        		for (Path p : ksp) {
            		outWriterAll.append(s +" "+target +" "+ p.toString()+"\n");
            	}
        		source_count+=1;
        		cur_k=k;
        		sum=0;
        		prev_sum=0;
        	}
        	else {
        		cur_k*=10;
//        		System.out.println(cur_k);
        		prev_sum=sum;
        		sum=0;
        	}
 
        }
        ksp=null;
        graph=null;
        outWriterAll.close();
        
        System.out.println(target);
        
        
    }
    
    public static void shortest_paths_by_length_on_unweighted(String graphFilename, String outFilename,
    		String source, String target, int k, int len_threshold) throws IOException
    {
        
    	/* Read graph from file */
        Graph graph = new Graph(graphFilename);

        Eppstein eppsteinAlgorithm = new Eppstein();
        
//      source is a string with multiple comma separated sources
        source = source.replace("[", "");
        source = source.replace("]", "");
        String[] sources = source.split("," ,-2);
//        System.out.println("Sources: "+ sources[0]);
        PrintWriter outWriterAll = new PrintWriter(new FileOutputStream(new File(outFilename), false));
    	List<Path> ksp;
    	
    	int source_count = 0;
    	int cur_k = k;
    	
		
        while (source_count<sources.length) {
    		HashMap<Integer, Integer> paths_per_len = new HashMap<Integer, Integer>();
        	String s = sources[source_count];
        	s = s.replace(" ","");
//        	System.out.println("s: "+ s);
        	ksp = eppsteinAlgorithm.ksp(graph, s, target, cur_k);
        	
        	//check if the fraction of contributions from shortest paths sum up to (1-tolerance). If not then run ksp with higher K.
        	int cur_path_len = 0;
        	for (Path p : ksp) {
        		cur_path_len = p.size();
        		if (paths_per_len.containsKey(cur_path_len)){
            		paths_per_len.put(cur_path_len,paths_per_len.get(cur_path_len)+1 );
        		}
        		else {
            		paths_per_len.put(cur_path_len, 1);
        		}
        	}
        	
        	if (cur_path_len >= (len_threshold+1)) { //This means the code exhausted over all the paths of length len_threshold and we can stop increasing k
        		
        		
//            	System.out.println("cur_k "+ cur_k);
//        		System.out.println(paths_per_len.toString());
        		outWriterAll.append(s +"\t"+target +"\t"+ paths_per_len.toString()+"\n");
//        		
//        		for (Path p : ksp) {
//            		outWriterAll.append(s +" "+target +" "+ p.toString()+"\n");
//            	}
        		source_count+=1;
        		cur_k=k;
        		}
        	else {
//        		System.out.println("cur_k "+ cur_k);
//            	System.out.println(paths_per_len.toString());
        		cur_k*=5;
            	

        	}
        	
        	
        }
        	
        	
        ksp=null;
        graph=null;
        outWriterAll.close();
        
//        System.out.println(target);
        
        
    }
    
    public static HashMap<Integer,Integer> parse_string_to_hashmap(HashMap<String, String> pathleninfo_per_source, String source) {
    	// Parse String to HasMap
		HashMap<Integer, Integer> n_paths_per_len = new HashMap<Integer,Integer>();
		String path_len_description = pathleninfo_per_source.get(source);
		path_len_description = path_len_description.replace("{", "");
		path_len_description = path_len_description.replace("}", "");
		path_len_description = path_len_description.replace(" ", "");
		
		String[] path_len_descriptions = path_len_description.split(",");
		  
		for (String path_desc : path_len_descriptions) {//path_desc = "3=33"
			n_paths_per_len.put(Integer.parseInt(path_desc.split("=")[0]), Integer.parseInt(path_desc.split("=")[1]));
		     }
		return n_paths_per_len;
      	
    }

    public static void shortest_paths_by_length_on_weighted(String graphFilename, String outFilename,
    		String source, String target, int k, int len_threshold, String len_threshold_filename) throws IOException
    {
    	
//    	Read how many paths are there for what length for a particular source-target pair
		HashMap<String, String> pathleninfo_per_source = new HashMap<String, String>();
    	try {
            BufferedReader in = new BufferedReader(new FileReader(len_threshold_filename));
            String line = in.readLine();
            while (line != null) {
            	String cur_source = line.split("\t")[0];
            	String cur_target = line.split("\t")[1];
            	pathleninfo_per_source.put(cur_source, line.split("\t")[2]); // path_len_description = "{3=33, 4=2681, 5=2286}"
                line = in.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    	
    	/* Read graph from file */
        Graph graph = new Graph(graphFilename);

        Eppstein eppsteinAlgorithm = new Eppstein();
        
//      source is a string with multiple comma separated sources
        source = source.replace("[", "");
        source = source.replace("]", "");
        String[] sources = source.split("," ,-2);
     
        /* Output the K shortest paths */
        PrintWriter outWriterAll = new PrintWriter(new FileOutputStream(new File(outFilename), false));
    	List<Path> ksp;
    	
    	int source_count = 0;
    	int cur_k = k;
    	
		
        while (source_count<sources.length) {
        	int flag = 0;
        	String s = sources[source_count];
        	s = s.replace(" ","");
//        	System.out.println("s: "+ s);
        	
    		HashMap<Integer, Integer> threshold_paths_per_len = parse_string_to_hashmap(pathleninfo_per_source, s);
    		HashMap<Integer, Integer> found_paths_per_len = new  HashMap<Integer, Integer>();
    		for (int i=0; i< len_threshold +2 ; i++ ) {
    			found_paths_per_len.put(i, 0);  			
    		}



        	ksp = eppsteinAlgorithm.ksp(graph, s, target, cur_k);
        	
        	//check if the fraction of contributions from shortest paths sum up to (1-tolerance). If not then run ksp with higher K.
        	for (Path p : ksp) {
//        		System.out.println(p.toString());
        		int cur_path_len = p.size();
        		if (found_paths_per_len.containsKey(cur_path_len)){
        			found_paths_per_len.put(cur_path_len,found_paths_per_len.get(cur_path_len)+1 );
        		}
        		else {
        			found_paths_per_len.put(cur_path_len, 1);
        		}
        	}
        	
//        	compare between threshold_paths_per_len and found_paths_per_len
//        	threshold_paths_per_len.forEach((key,value)->System.out.println(key + " = " + value) );
        	
        	for (Map.Entry<Integer, Integer> set : threshold_paths_per_len.entrySet()) {
        		Integer  th_path_len = set.getKey();
        		Integer th_n_paths = set.getValue();
        		Integer found_paths = found_paths_per_len.get(th_path_len);
        		
        		if ((th_path_len<=len_threshold)) {
        			if ((found_paths <th_n_paths)) {
        				flag = 1; //this means we need to run for larger k
            			break;
        			}
            	}
            }
        	
        	if (flag==1) {
//        		System.out.println("cur_k "+ cur_k);
//        		System.out.println("Path needed: " + threshold_paths_per_len.toString());
//            	System.out.println("Path found: " + found_paths_per_len.toString());

        		cur_k*=5;
        	}
        	
        	else {
//        		System.out.println("cur_k "+ cur_k);
//            	System.out.println("Path needed: " + threshold_paths_per_len.toString());
//            	System.out.println("Path found: " + found_paths_per_len.toString());

        		for (Path p : ksp) {
            		outWriterAll.append(s +" "+target +" "+ p.toString()+"\n");
            	}
        		source_count+=1;
        		cur_k=k;
        		
        	}
        }
        	
        	
        ksp=null;
        graph=null;
        outWriterAll.close();
//        System.out.println(target);
        
        
    }
    
    
}
