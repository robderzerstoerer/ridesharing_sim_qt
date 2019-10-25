#include <cstddef>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <queue>
#include <map>
#include <list>
#include <set>
#include <tuple>
#include <vector>
#include <cmath>
#include <random>

#include <cassert>

#include "ridesharing_sim.h"

#include "network.hpp"
#include "delaunay_network.h"

#ifndef _INTEGER_TYPES
	#define ULL uint64_t
	#define LL int64_t
	#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
	#define MACRO_EPSILON 0.000000000001
	#define _EPSILON
#endif

constexpr double pi = 3.14159265358979323846;

//ULL and double: 64-bit variables for everything

//Simulation of taxi system with arbitrary networks (needs to be strongly connected) and taxis with different service types
//everything is independent of the destination of a request (constant maximal waiting time for all requests, indiscriminate service of all requests)
int main(int argc, char* argv[])
{
//	ULL number_of_buses = 25;
	ULL number_of_nodes = 100;
	ULL number_of_poisson_random_edges = 200;
//
//	double normalized_request_rate = 1;
//	std::string topology = "ring";

//	std::vector<std::string> topology_list = {"ring", "torus"};
//	std::vector<std::string> topology_list = {"two_nodes"};
//	std::vector<std::string> topology_list = {"ring"};
//	std::vector<std::string> topology_list = {"directed_ring"};
//	std::vector<std::string> topology_list = {"torus"};
//	std::vector<std::string> topology_list = {"double_ring"};
//	std::vector<std::string> topology_list = {"ladder"};
//	std::vector<std::string> topology_list = {"3ladder"};
//	std::vector<std::string> topology_list = {"4ladder"};
//	std::vector<std::string> topology_list = {"star"};	//needs much linger equilibration time
//	std::vector<std::string> topology_list = {"3cayley_tree"};
//	std::vector<std::string> topology_list = {"simplified_city"};
//	std::vector<std::string> topology_list = {"complete_graph"};
//	std::vector<std::string> topology_list = {"poisson_random"};
	std::vector<std::string> topology_list = {"delaunay_random_torus"};
//	std::vector<std::string> topology_list = {"delaunay_random_grid"};
//	std::vector<std::string> topology_list = {"gabriel_random_torus"};
//	std::vector<std::string> topology_list = {"gabriel_random_grid"};
//	std::vector<ULL> number_of_buses_list = {10};
//	std::vector<ULL> number_of_buses_list = {0,1,2,3,4,5,6,7,8,9,10,100,1000};
	std::vector<ULL> number_of_buses_list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000, 5000};//,6000,7000,8000,10000};
//	std::vector<ULL> number_of_buses_list = {600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000, 5000};//,6000,7000,8000,10000};
//	std::vector<ULL> number_of_buses_list = {7500, 10000, 12500, 15000, 20000, 25000};
//	std::vector<ULL> number_of_buses_list = {2, 10, 20, 50, 100, 200, 500};
//	std::vector<double> normalized_request_rate_list = {0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0};
//	std::vector<double> normalized_request_rate_list = {0.5,1,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10};
//	std::vector<double> normalized_request_rate_list = {0.5,1,2,5,10};
	std::vector<double> normalized_request_rate_list = {15};

//	std::vector<ULL> number_of_buses_list = {2, 10, 100, 500};
//	std::vector<double> normalized_request_rate_list = {0.5,2};

	std::stringstream filename("");
//	filename << "test.dat";
//	filename << "min_arrival_time__uniform_requests__delaunay_random_all_dx_025.dat";
	filename << "min_arrival_time__uniform_requests__delaunay_N_100__x_15__corrected_assignment_test.dat";
	std::ofstream out(filename.str().c_str());


	for(std::string topology : topology_list)
	{
//		std::stringstream filename("");
//		filename << "./RESULTS__min_arrival_time__uniform_requests/min_arrival_time__uniform_requests__" << topology << "_N_" << number_of_nodes << ".dat";
//		std::ofstream out(filename.str().c_str());

		for(ULL number_of_buses : number_of_buses_list)
		{
			for(double normalized_request_rate : normalized_request_rate_list)
			{

				ridesharing_sim sim(number_of_nodes, number_of_buses, 0);

				for(unsigned int i = 0; i < number_of_nodes; ++i)
				{
					if(topology == "two_nodes" )
					{
						assert(number_of_nodes == 2);

						if(i == 0)
						{
							sim.network.add_link(0, 1, 1);
							sim.network.add_link(1, 0, 1);
						}
					}

					if(topology == "line" || topology == "ring")
					{
						//line
						if(i < number_of_nodes - 1)
						{
							sim.network.add_link(i  , i+1, 1);
							sim.network.add_link(i+1, i  , 1);
						}

						//ring
						if(topology == "ring")
						{
							if(i == number_of_nodes-1)
							{
								sim.network.add_link(0  , number_of_nodes - 1, 1);
								sim.network.add_link(number_of_nodes - 1, 0  , 1);
							}
						}
					}

					if(topology == "directed_ring")
					{
						//line
						if(i < number_of_nodes - 1)
						{
							sim.network.add_link(i  , i+1, 1);
						}

						if(i == number_of_nodes-1)
						{
							sim.network.add_link(number_of_nodes - 1, 0  , 1);
						}
					}

					if(topology == "double_ring")
					{
						if(i < number_of_nodes/2 - 1)
						{
							sim.network.add_link(i  , i+1, 1);
							sim.network.add_link(i+1, i  , 1);
						}else if(i > number_of_nodes/2 - 1 && i < number_of_nodes - 1)
						{
							sim.network.add_link(i  , i+1, 1);
							sim.network.add_link(i+1, i  , 1);
						}

						if(i == number_of_nodes/2 - 1)
						{
							sim.network.add_link(0 , number_of_nodes/2 - 1, 1);
							sim.network.add_link(number_of_nodes/2 - 1, 0 , 1);
						}

						if(i == number_of_nodes - 1)
						{
							sim.network.add_link(number_of_nodes/2  , number_of_nodes - 1, 1);
							sim.network.add_link(number_of_nodes - 1, number_of_nodes/2  , 1);
						}

						if(i == 0)
						{
							sim.network.add_link(number_of_nodes/2  , 0, 1);
							sim.network.add_link(0, number_of_nodes/2  , 1);
						}
					}

					if(topology == "ladder")
					{
						if(i < number_of_nodes/2 - 1)
						{
							sim.network.add_link(i  , i+1, 1);
							sim.network.add_link(i+1, i  , 1);
						}else if(i > number_of_nodes/2 - 1 && i < number_of_nodes - 1)
						{
							sim.network.add_link(i  , i+1, 1);
							sim.network.add_link(i+1, i  , 1);
						}

						if(i == number_of_nodes/2 - 1)
						{
							sim.network.add_link(0 , number_of_nodes/2 - 1, 1);
							sim.network.add_link(number_of_nodes/2 - 1, 0 , 1);
						}

						if(i == number_of_nodes - 1)
						{
							sim.network.add_link(number_of_nodes/2  , number_of_nodes - 1, 1);
							sim.network.add_link(number_of_nodes - 1, number_of_nodes/2  , 1);
						}

						if(i < number_of_nodes/2)
						{
							sim.network.add_link(i + number_of_nodes/2  , i, 1);
							sim.network.add_link(i, i + number_of_nodes/2  , 1);
						}
					}

					if(topology == "3ladder")
					{
						unsigned int k = 3;
						for(unsigned int m = 0; m < k; ++m)
						{
							if(i >= m * number_of_nodes/k && i < (m+1)*number_of_nodes/k - 1)
							{
								sim.network.add_link(i  , i+1, 1);
								sim.network.add_link(i+1, i  , 1);

								sim.network.add_link(i  , (i+number_of_nodes/k)%number_of_nodes, 1);
								sim.network.add_link((i+number_of_nodes/k)%number_of_nodes, i  , 1);
							}

							if(i+1 == (m+1)*number_of_nodes/k )
							{
								sim.network.add_link(m*number_of_nodes/k , (m+1)*number_of_nodes/k - 1, 1);
								sim.network.add_link((m+1)*number_of_nodes/k - 1, m*number_of_nodes/k , 1);

								sim.network.add_link(i  , (i+number_of_nodes/k)%number_of_nodes, 1);
								sim.network.add_link((i+number_of_nodes/k)%number_of_nodes, i  , 1);
							}
						}
					}

					if(topology == "4ladder")
					{
						unsigned int k = 4;
						for(unsigned int m = 0; m < k; ++m)
						{
							if(i >= m * number_of_nodes/k && i < (m+1)*number_of_nodes/k - 1)
							{
								sim.network.add_link(i  , i+1, 1);
								sim.network.add_link(i+1, i  , 1);

								sim.network.add_link(i  , (i+number_of_nodes/k)%number_of_nodes, 1);
								sim.network.add_link((i+number_of_nodes/k)%number_of_nodes, i  , 1);
							}

							if(i+1 == (m+1)*number_of_nodes/k )
							{
								sim.network.add_link(m*number_of_nodes/k , (m+1)*number_of_nodes/k - 1, 1);
								sim.network.add_link((m+1)*number_of_nodes/k - 1, m*number_of_nodes/k , 1);

								sim.network.add_link(i  , (i+number_of_nodes/k)%number_of_nodes, 1);
								sim.network.add_link((i+number_of_nodes/k)%number_of_nodes, i  , 1);
							}
						}
					}

//					if(topology == "star")
//					{
//						//star
//						if(i < number_of_nodes-1)
//						{
//							sim.network.add_link(i  , number_of_nodes - 1, 1);
//							sim.network.add_link(number_of_nodes - 1, i  , 1);
//						}
//					}
//
					if(topology == "star")
					{
						//star
						if(i > 0)
						{
							sim.network.add_link(i  , 0, 1);
							sim.network.add_link(0, i  , 1);
						}
					}

					if(topology == "grid" || topology == "torus")
					{
						//grid
						unsigned int L = sqrt(number_of_nodes);
						if(i%L != L-1)
						{
							sim.network.add_link(i  , i+1, 1);
							sim.network.add_link(i+1, i  , 1);
						}
						if(i < L*(L-1))
						{
							sim.network.add_link(i  , i+L, 1);
							sim.network.add_link(i+L, i  , 1);
						}
						if(topology == "torus")
						{
							//torus
							if(i%L == L-1)
							{
								sim.network.add_link(i    , i+1-L, 1);
								sim.network.add_link(i+1-L, i    , 1);
							}
							if(i >= L*(L-1))
							{
								sim.network.add_link(i  , (i+L)%number_of_nodes, 1);
								sim.network.add_link((i+L)%number_of_nodes, i  , 1);
							}
						}
					}

					if(topology == "simplified_city")
					{
						//simplified city (spiderweb)
						unsigned int arms = sqrt(number_of_nodes);
						unsigned int nodes_per_arm = sqrt(number_of_nodes);
						unsigned int nodes_between_circles = nodes_per_arm/floor(sqrt(nodes_per_arm));
						//arms
						if( i%nodes_per_arm != 0)
						{
							sim.network.add_link(i  , i-1, 1);
							sim.network.add_link(i-1, i  , 1);
						}

						//rings
						if( (i%nodes_per_arm)%nodes_between_circles == 0 )
						{
							sim.network.add_link(i  , (i+nodes_per_arm)%number_of_nodes, 1);
							sim.network.add_link((i+nodes_per_arm)%number_of_nodes, i  , 1);
						}
					}

					if(topology == "complete_graph")
					{
						for(unsigned int j = i+1; j < number_of_nodes; ++j )
						{
							sim.network.add_link(i  , j, 1);
							sim.network.add_link(j, i  , 1);
						}
					}

					if(topology == "3cayley_tree")
					{
						for(unsigned int i = 0; i < number_of_nodes; ++i)
						{
							if(i > 0)
							{
								if(i < 4)
								{
									sim.network.add_link(i , 0, 1);
									sim.network.add_link(0 , i, 1);
								}else{
									sim.network.add_link(i, i/2-1, 1);
									sim.network.add_link(i/2-1, i, 1);
								}
							}
						}
					}

					if(topology == "poisson_random")
					{
                        poisson_random_network<std::mt19937_64> network_edgelist( sim.random_generator, number_of_nodes, number_of_poisson_random_edges );

						std::vector< std::vector< long int > > temp_dist(number_of_nodes, std::vector<long int >(number_of_nodes, 100000));

                        while( !distances(network_edgelist.network_edgelist, number_of_nodes, temp_dist) )
						{
							network_edgelist.random_network_links(number_of_nodes, number_of_poisson_random_edges);
						}

                        for(auto n : network_edgelist.edges())
							for(auto e : n.second)
								sim.network.add_link( n.first, e, 1 );

						break;
					}

					if(topology == "delaunay_random_torus")
					{
						delaunay_network<std::mt19937_64> delaunay_net(sim.random_generator, number_of_nodes, true);
						for(auto e : delaunay_net.edges())
						{
							sim.network.add_link( e.first, e.second, delaunay_net.distance(e.first, e.second) );
							sim.network.add_link( e.second, e.first, delaunay_net.distance(e.second, e.first) );
						}
						break;

//						//test torus topology
//						std::ofstream testout("delaunay_pos.dat");
//						for(unsigned int i = 0; i < number_of_nodes; ++i)
//							testout << i << '\t' << delaunay_net.position(i).x() << '\t' << delaunay_net.position(i).y() << std::endl;
//						testout.close();
//
//						testout.open("delaunay_links.dat");
//						for(auto e : delaunay_net.edges())
//						{
////							if(delaunay_net.distance(e.first, e.second) < 0.2)
//								testout << e.first << '\t' << e.second << std::endl;
//						}
//						testout.close();
//
//						return(1);
					}

					if(topology == "gabriel_random_torus")
					{
						delaunay_network<std::mt19937_64> delaunay_net(sim.random_generator, number_of_nodes, true);
						delaunay_net.random_gabriel_network();

						for(auto e : delaunay_net.edges())
						{
							sim.network.add_link( e.first, e.second, delaunay_net.distance(e.first, e.second) );
							sim.network.add_link( e.second, e.first, delaunay_net.distance(e.second, e.first) );
						}

						break;
					}
				}

			//	ULL temp_from, temp_to;
			//	double temp_dist;
			//	std::ifstream in("test_edgelist_Goettingen_converted.dat");
			//	while( !in.eof() )
			//	{
			//		in >> temp_from >> temp_to >> temp_dist;
			////		edgelist[temp_from].insert( std::make_pair(temp_to, temp_dist) );
			//		edgelist[temp_from].insert( std::make_pair(temp_to, 1) );
			//	}


				//output network topology
//				std::stringstream ss;
//				ss << topology << "_N_" << number_of_nodes << "_edgelist.dat";
//				std::ofstream testout(ss.str().c_str());
//				for(auto i : sim.network.edgelist)
//					for(auto e : i.second)
//						testout << i.first << '\t' << e.first << '\t' << e.second << std::endl;
//				testout.close();
//				return(1);



				sim.network.create_distances();

				for(unsigned int i = 0; i < number_of_nodes; ++i)
					for(unsigned int j = 0; j < number_of_nodes; ++j)
						assert(sim.network.get_network_distance(i,j) < 10 * number_of_nodes);

				//set probabilities for requests
				std::vector<double> origin_probabilities(number_of_nodes, 1.0 / number_of_nodes);
				std::vector<double> destination_probabilities(number_of_nodes, 1.0 / number_of_nodes);

				std::uniform_real_distribution<double> uniform01 (0, 1);				//generic uniform distribution
				std::uniform_int_distribution<ULL> random_node(0, number_of_nodes-1);	//uniform distribution over all nodes

//				//random request
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					origin_probabilities[i] = uniform01(sim.random_generator);
//					destination_probabilities[i] = origin_probabilities[i];
//				}

//				//random request asymmetric
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					origin_probabilities[i] = uniform01(sim.random_generator);
//					destination_probabilities[i] = uniform01(sim.random_generator);
//				}

				//gaussian requests
//				ULL origin_center = 0;
//				ULL destination_center = origin_center;
//				double max_dist = 0;
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					if( sim.network.get_network_distance(i,origin_center) > max_dist )
//						max_dist = sim.network.get_network_distance(i,origin_center);
//				}
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					origin_probabilities[i] = exp( - 2 * sim.network.get_network_distance(i,origin_center)*sim.network.get_network_distance(i,origin_center) / ( max_dist*max_dist ) );
//					destination_probabilities[i] = origin_probabilities[i];
//				}

//				//gaussian requests asymmetric
//				ULL origin_center = 0;
////				ULL destination_center = 55;	//specifically for torus N=100
//				ULL destination_center = 61;	//specifically for torus N=121
//				double max_dist_origin = 0;
//				double max_dist_dest = 0;
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					if( sim.network.get_network_distance(i,origin_center) > max_dist_origin )
//						max_dist_origin = sim.network.get_network_distance(i,origin_center);
//					if( sim.network.get_network_distance(i,destination_center) > max_dist_dest )
//						max_dist_dest = sim.network.get_network_distance(i,destination_center);
//				}
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					origin_probabilities[i] = exp( - 2 * sim.network.get_network_distance(i,origin_center)*sim.network.get_network_distance(i,origin_center) / ( max_dist_origin*max_dist_origin ) );
//					destination_probabilities[i] = exp( - 2 * sim.network.get_network_distance(i,destination_center)*sim.network.get_network_distance(i,destination_center) / ( max_dist_dest*max_dist_dest ) );
//				}

				//gaussian requests strongly asymmetric
//				ULL origin_center = 0;
////				ULL destination_center = 55;	//specifically for torus N=100
//				ULL destination_center = 61;	//specifically for torus N=121
//				double max_dist_origin = 0;
//				double max_dist_dest = 0;
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					if( sim.network.get_network_distance(i,origin_center) > max_dist_origin )
//						max_dist_origin = sim.network.get_network_distance(i,origin_center);
//					if( sim.network.get_network_distance(i,destination_center) > max_dist_dest )
//						max_dist_dest = sim.network.get_network_distance(i,destination_center);
//				}
//				for( ULL i = 0; i < number_of_nodes; ++i )
//				{
//					origin_probabilities[i] = exp( - 8 * sim.network.get_network_distance(i,origin_center)*sim.network.get_network_distance(i,origin_center) / ( max_dist_origin*max_dist_origin ) );
//					destination_probabilities[i] = exp( - 8 * sim.network.get_network_distance(i,destination_center)*sim.network.get_network_distance(i,destination_center) / ( max_dist_dest*max_dist_dest ) );
//				}


				//reset the
//				sim.network.set_origin_probabilities();	//default: uniform
//				sim.network.set_destination_probabilities();	//default: uniform
				sim.network.set_origin_probabilities(origin_probabilities);
				sim.network.set_destination_probabilities(destination_probabilities);


				sim.network.recalc_mean_distances();

			//	std::cout << sim.network.get_mean_pickup_distance() << '\t' << sim.network.get_mean_dropoff_distance() << std::endl;

				for(ULL b = 0; b < number_of_buses; ++b)
				{
			//		std::cout << b << '\t' << sim.transporter_list[b].get_index() << '\t' << sim.transporter_list[b].get_current_location() << std::endl;

					sim.transporter_list[b].reset( b, sim.network.generate_request().second, 0 );

			//		std::cout << b << '\t' << sim.transporter_list[b].get_index() << '\t' << sim.transporter_list[b].get_current_location() << std::endl;
				}

			//	sim.set_normalized_request_rate(3);
			//	sim.run_sim_time(1000);
			//	sim.enable_measurements(10);
			//	sim.enable_timeseries_output(0.1, "time_series_test.dat");
			//	sim.run_sim_requests(100*number_of_buses);
			//	sim.run_sim_time(10000);

				std::cout << "running simulation: " << topology << ", B = " << number_of_buses << ", x = " << normalized_request_rate <<std::endl;

				sim.set_normalized_request_rate(normalized_request_rate);
//				sim.set_diaptcher({"min_arrival","max_pickup","max_occupancy"});

				//use this to ensure uniform distribution of transporters for the two_nodes topology
//				std::list< std::pair< double, std::pair<ULL,ULL> > > request_list;
//				for(double t = 0; t < 10.0; t += 2.0 / number_of_buses)
//				{
//					request_list.push_back( std::make_pair( t, std::make_pair(0,1) ) );
//					request_list.push_back( std::make_pair( t, std::make_pair(1,0) ) );
//				}


//				//use this to ensure uniform distribution of transporters for the ring topology
//				std::list< std::pair< double, std::pair<ULL,ULL> > > request_list;
//				for(double t = 0; t < 10.0; t += 2.0 * number_of_nodes / (1.0 * number_of_buses) )
//				{
//					for(ULL i = 0; i < number_of_nodes; ++i)
//					{
//						request_list.push_back( std::make_pair( t, std::make_pair(i, (i+1) % number_of_nodes) ) );
//						request_list.push_back( std::make_pair( t, std::make_pair((i+1) % number_of_nodes, i) ) );
//					}
//				}

//				//use this to ensure uniform distribution of transporters for the four node star topology
//				std::list< std::pair< double, std::pair<ULL,ULL> > > request_list;
//				for(double t = 0; t < 6.0 * (number_of_nodes-1.0) / (long double) number_of_buses; )//t += 2.0 * (number_of_nodes-1.0) / (long double) number_of_buses)
//					//t from 0 to 6 trips along edge, distance between buses:
//				{
//					for(ULL i = 1; i < number_of_nodes; ++i)
//					{
//						request_list.push_back( std::make_pair( t, std::make_pair(0,i) ) );
//						request_list.push_back( std::make_pair( t, std::make_pair(i,0) ) );
//
//						t += 2.0 / number_of_buses;
//					}
//				}

////				//use this to ensure uniform distribution of transporters for the complete graph topology
//				std::list< std::pair< double, std::pair<ULL,ULL> > > request_list;
//				for(double t = 0; t < 5.0; t += (1.0 * number_of_nodes*(number_of_nodes-1)) / number_of_buses)
//				{
//					for(unsigned int i = 0; i < number_of_nodes; ++i)
//						for(unsigned int j = 0; j < number_of_nodes; ++j)
//							if(i != j)
//								request_list.push_back( std::make_pair( t, std::make_pair(i,j) ) );
//				}
//
//				use this to ensure uniform distribution of transporters for a torus
				std::list< std::pair< double, std::pair<ULL,ULL> > > request_list;
				if(number_of_buses > 0.25 * number_of_nodes)
				{
					for(double t = 0; t < 1.0; t += (4.0 * number_of_nodes) / (number_of_buses) )
					{
						unsigned int L = sqrt(number_of_nodes);

						for(unsigned int i = 0; i < number_of_nodes; ++i)
						{
							if(i % L == 0)
								request_list.push_back( std::make_pair( t, std::make_pair(i, i+L-1) ) );
							else
								request_list.push_back( std::make_pair( t, std::make_pair(i, i-1) ) );

							if(i % L == L-1)
								request_list.push_back( std::make_pair( t, std::make_pair(i, i+1-L) ) );
							else
								request_list.push_back( std::make_pair( t, std::make_pair(i, i+1) ) );

							if(i + L >= number_of_nodes )
								request_list.push_back( std::make_pair( t, std::make_pair(i, i+L-number_of_nodes) ) );
							else
								request_list.push_back( std::make_pair( t, std::make_pair(i, i+L) ) );

							if(i < L)
								request_list.push_back( std::make_pair( t, std::make_pair(i, i+number_of_nodes-L) ) );
							else
								request_list.push_back( std::make_pair( t, std::make_pair(i, i-L) ) );
						}
					}
				}

//				std::cout << request_list.size() << std::endl;

				sim.run_sim_request_list( request_list );

				sim.set_normalized_request_rate(normalized_request_rate);
				sim.run_sim_requests( std::max((ULL)10000, 100*number_of_buses) );
				sim.enable_measurements( std::max(1.123456789,1.123456789/sim.request_rate) );
				sim.run_sim_requests( std::max((ULL)100000, 1000*number_of_buses) );
//				sim.run_sim_requests( std::max((ULL)100000, 100*number_of_buses) );
//				sim.enable_measurements( ((number_of_buses <= 10)?(1.0):(0.1))*pi/normalized_request_rate );	//pi/(10 * normalized request rate) for 100 or 1000 buses due to larger actual request rates
//				sim.run_sim_requests( 1e9 );

//				out << std::setprecision(10);
				sim.print_params(out);
				sim.print_measurements(out);

//				std::cin.ignore();
			}
		}
//		out.close();
	}

	out.close();

	return(0);
}

