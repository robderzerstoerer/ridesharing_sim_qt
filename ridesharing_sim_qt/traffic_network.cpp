#include "traffic_network.h"
#include "network.h"
#include <functional>

traffic_network::traffic_network()
{
	random_generator.seed(1);
}

traffic_network::traffic_network(ULL param_N)
{
	random_generator.seed(1);
	init(param_N, std::vector< std::tuple<ULL, ULL, double> >());	//call init with an empty vector of links
}

traffic_network::traffic_network(ULL param_N, std::vector< std::tuple<ULL, ULL, double> > param_links)
{
	random_generator.seed(1);
	init(param_N, param_links);
}

//initialize the network and the relevant data structures
void traffic_network::init(ULL param_N, std::vector< std::tuple<ULL, ULL, double> > param_links)
{
	number_of_nodes = param_N;
	uniform01 = std::uniform_real_distribution<double>(0, 1);

	potential_route_nodes.reserve(number_of_nodes);

	network_distances.resize(number_of_nodes);
	for (ULL i = 0; i < number_of_nodes; ++i)
	{
		network_distances[i] = std::vector<double>(number_of_nodes, 1e10);
	}

	//add all links in the list
	for (auto& e : param_links)
		add_link(std::get<0>(e), std::get<1>(e), std::get<2>(e));

	//initialize default value for probabilities: uniform distribution
	origin_probabilities = std::vector<double>(number_of_nodes, 1.0 / number_of_nodes);
	random_origin = std::discrete_distribution<ULL>(origin_probabilities.begin(), origin_probabilities.end());
	destination_probabilities = std::vector<double>(number_of_nodes, 1.0 / number_of_nodes);
	random_destination = std::discrete_distribution<ULL>(destination_probabilities.begin(), destination_probabilities.end());
	recalc_mean_distances();
}

void traffic_network::init(std::string topology, ULL param_N)
{
	random_generator.seed(1);
	init(param_N, std::vector< std::tuple<ULL, ULL, double> >());
	//create topologies
	//add links for each node in the network
	for (unsigned int i = 0; i < number_of_nodes; ++i)
	{
		if (topology == "two_nodes")
		{
			assert(param_N == 2);

			if (i == 0)
			{
				add_link(0, 1, 1);
				add_link(1, 0, 1);
			}
		}

		if (topology == "ring")
		{
			//create a path-graph (a line)
			if (i < number_of_nodes - 1)
			{
				add_link(i, i + 1, 1);
				add_link(i + 1, i, 1);
			}

			//close the ring
			if (i == number_of_nodes - 1)
			{
				add_link(0, number_of_nodes - 1, 1);
				add_link(number_of_nodes - 1, 0, 1);
			}

		}

		if (topology == "directed ring")
		{
			//create a path-graph (a line)
			if (i < number_of_nodes - 1)
			{
				add_link(i, i + 1, 1);
			}

			//close the ring
			if (i == number_of_nodes - 1)
			{
				add_link(number_of_nodes - 1, 0, 1);
			}

		}

		if (topology == "double_ring")
		{
			if (i < number_of_nodes / 2 - 1)
			{
				add_link(i, i + 1, 1);
				add_link(i + 1, i, 1);
			}
			else if (i > number_of_nodes / 2 - 1 && i < number_of_nodes - 1)
			{
				add_link(i, i + 1, 1);
				add_link(i + 1, i, 1);
			}

			if (i == number_of_nodes / 2 - 1)
			{
				add_link(0, number_of_nodes / 2 - 1, 1);
				add_link(number_of_nodes / 2 - 1, 0, 1);
			}

			if (i == number_of_nodes - 1)
			{
				add_link(number_of_nodes / 2, number_of_nodes - 1, 1);
				add_link(number_of_nodes - 1, number_of_nodes / 2, 1);
			}

			if (i == 0)
			{
				add_link(number_of_nodes / 2, 0, 1);
				add_link(0, number_of_nodes / 2, 1);
			}
		}

		if (topology == "ladder")
		{
			if (i < number_of_nodes / 2 - 1)
			{
				add_link(i, i + 1, 1);
				add_link(i + 1, i, 1);
			}
			else if (i > number_of_nodes / 2 - 1 && i < number_of_nodes - 1)
			{
				add_link(i, i + 1, 1);
				add_link(i + 1, i, 1);
			}

			if (i == number_of_nodes / 2 - 1)
			{
				add_link(0, number_of_nodes / 2 - 1, 1);
				add_link(number_of_nodes / 2 - 1, 0, 1);
			}

			if (i == number_of_nodes - 1)
			{
				add_link(number_of_nodes / 2, number_of_nodes - 1, 1);
				add_link(number_of_nodes - 1, number_of_nodes / 2, 1);
			}

			if (i < number_of_nodes / 2)
			{
				add_link(i + number_of_nodes / 2, i, 1);
				add_link(i, i + number_of_nodes / 2, 1);
			}
		}

		if (topology == "3ladder")
		{
			unsigned int k = 3;
			for (unsigned int m = 0; m < k; ++m)
			{
				if (i >= m * number_of_nodes / k && i < (m + 1) * number_of_nodes / k - 1)
				{
					add_link(i, i + 1, 1);
					add_link(i + 1, i, 1);

					add_link(i, (i + number_of_nodes / k) % number_of_nodes, 1);
					add_link((i + number_of_nodes / k) % number_of_nodes, i, 1);
				}

				if (i + 1 == (m + 1) * number_of_nodes / k)
				{
					add_link(m * number_of_nodes / k, (m + 1) * number_of_nodes / k - 1, 1);
					add_link((m + 1) * number_of_nodes / k - 1, m * number_of_nodes / k, 1);

					add_link(i, (i + number_of_nodes / k) % number_of_nodes, 1);
					add_link((i + number_of_nodes / k) % number_of_nodes, i, 1);
				}
			}
		}

		if (topology == "4ladder")
		{
			unsigned int k = 4;
			for (unsigned int m = 0; m < k; ++m)
			{
				if (i >= m * number_of_nodes / k && i < (m + 1) * number_of_nodes / k - 1)
				{
					add_link(i, i + 1, 1);
					add_link(i + 1, i, 1);

					add_link(i, (i + number_of_nodes / k) % number_of_nodes, 1);
					add_link((i + number_of_nodes / k) % number_of_nodes, i, 1);
				}

				if (i + 1 == (m + 1) * number_of_nodes / k)
				{
					add_link(m * number_of_nodes / k, (m + 1) * number_of_nodes / k - 1, 1);
					add_link((m + 1) * number_of_nodes / k - 1, m * number_of_nodes / k, 1);

					add_link(i, (i + number_of_nodes / k) % number_of_nodes, 1);
					add_link((i + number_of_nodes / k) % number_of_nodes, i, 1);
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
		if (topology == "star")
		{
			//star
			if (i > 0)
			{
				add_link(i, 0, 1);
				add_link(0, i, 1);
			}
		}

		if (topology == "grid" || topology == "torus")
		{
			//grid
			unsigned int L = sqrt(number_of_nodes);
			if (i % L != L - 1)
			{
				add_link(i, i + 1, 1);
				add_link(i + 1, i, 1);
			}
			if (i < L * (L - 1))
			{
				add_link(i, i + L, 1);
				add_link(i + L, i, 1);
			}
			if (topology == "torus")
			{
				//torus
				if (i % L == L - 1)
				{
					add_link(i, i + 1 - L, 1);
					add_link(i + 1 - L, i, 1);
				}
				if (i >= L * (L - 1))
				{
					add_link(i, (i + L) % number_of_nodes, 1);
					add_link((i + L) % number_of_nodes, i, 1);
				}
			}
		}

		if (topology == "simplified_city")
		{
			//simplified city (spiderweb)
			unsigned int arms = sqrt(number_of_nodes);
			unsigned int nodes_per_arm = sqrt(number_of_nodes);
			unsigned int nodes_between_circles = nodes_per_arm / floor(sqrt(nodes_per_arm));
			//arms
			if (i % nodes_per_arm != 0)
			{
				add_link(i, i - 1, 1);
				add_link(i - 1, i, 1);
			}

			//rings
			if ((i % nodes_per_arm) % nodes_between_circles == 0)
			{
				add_link(i, (i + nodes_per_arm) % number_of_nodes, 1);
				add_link((i + nodes_per_arm) % number_of_nodes, i, 1);
			}
		}

		if (topology == "complete_graph")
		{
			for (unsigned int j = i + 1; j < number_of_nodes; ++j)
			{
				add_link(i, j, 1);
				add_link(j, i, 1);
			}
		}

		if (topology == "3cayley_tree")
		{
			for (unsigned int i = 0; i < number_of_nodes; ++i)
			{
				if (i > 0)
				{
					if (i < 4)
					{
						add_link(i, 0, 1);
						add_link(0, i, 1);
					}
					else {
						add_link(i, i / 2 - 1, 1);
						add_link(i / 2 - 1, i, 1);
					}
				}
			}
		}

		if (topology == "poisson_random")
		{
			poisson_random_network<std::mt19937_64> network_edgelist(random_generator, number_of_nodes, number_of_poisson_random_edges);

			std::vector< std::vector< long int > > temp_dist(number_of_nodes, std::vector<long int >(number_of_nodes, 100000));

			while (!distances(network_edgelist.network_edgelist, number_of_nodes, temp_dist))
			{
				network_edgelist.random_network_links(number_of_nodes, number_of_poisson_random_edges);
			}

			for (auto n : network_edgelist.edges())
				for (auto e : n.second)
					add_link(n.first, e, 1);

			break;
		}

		if (topology == "delaunay_random_torus")
		{
			delaunay_network<std::mt19937_64> delaunay_net(random_generator, number_of_nodes, true);
			for (auto e : delaunay_net.edges())
			{
				add_link(e.first, e.second, delaunay_net.distance(e.first, e.second));
				add_link(e.second, e.first, delaunay_net.distance(e.second, e.first));
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

		if (topology == "gabriel_random_torus")
		{
			delaunay_network<std::mt19937_64> delaunay_net(random_generator, number_of_nodes, true);
			delaunay_net.random_gabriel_network();

			for (auto e : delaunay_net.edges())
			{
				add_link(e.first, e.second, delaunay_net.distance(e.first, e.second));
				add_link(e.second, e.first, delaunay_net.distance(e.second, e.first));
			}

			break;
		}
	}

	//pre-calculate distance matrix (to find shortest paths later)
	create_distances();

	//set probability distribution for origin and destination of requests
	//requests are uncorrelated from one random node to another random node
	//note: this currently allows requests from a node to itself (this should probably be changed due to mathematical problems when trying to calculate with this in extreme cases)
	set_origin_probabilities();			//default: uniform
	set_destination_probabilities();	//default: uniform
													//recompute mean distance with respect to request distribution
	recalc_mean_distances();
}

//destructor, do all the cleanup
traffic_network::~traffic_network()
{
	potential_route_nodes.clear();
	route.clear();
	origin_probabilities.clear();
	destination_probabilities.clear();

	for (ULL i = 0; i < number_of_nodes; ++i)
	{
		network_distances[i].clear();
	}
	network_distances.clear();

	for (auto& e : edgelist)
		e.second.clear();
	edgelist.clear();
}

void traffic_network::add_link(ULL from, ULL to, double dist)
{
	//assert parameter bounds
	assert(from < number_of_nodes);
	assert(to < number_of_nodes);
	assert(dist >= 0);

	//links are directed!
	edgelist[from].insert(std::make_pair(to, dist));
}

//calculate all shortest path distances by breadth first search (direct implementation of Dijkstra)
void traffic_network::create_distances()
{
	//reset the shape for the distance matrix, with 1e10 distances between all nodes (no distances found yet)
	network_distances.resize(number_of_nodes);
	for (ULL i = 0; i < number_of_nodes; ++i)
	{
		network_distances[i] = std::vector<double>(number_of_nodes, 1e10);
	}

	//then fill the distance matrix by finding all shortest paths (works for positive weighted graphs)
	std::priority_queue< std::pair<double, ULL>, std::vector< std::pair<double, ULL> >, std::greater< std::pair<double, ULL> > > next;
	std::pair<double, ULL> current;

	for (ULL i = 0; i < number_of_nodes; ++i)
	{
		network_distances[i][i] = 0;
		next.push(std::make_pair(network_distances[i][i], i));

		while (!next.empty())
		{
			current = next.top();
			next.pop();

			if (network_distances[i][current.second] == current.first)
			{
				for (std::pair<ULL, double> e : edgelist[current.second])
				{
					if (network_distances[i][e.first] > network_distances[i][current.second] + e.second)
					{
						network_distances[i][e.first] = network_distances[i][current.second] + e.second;
						next.push(std::make_pair(network_distances[i][e.first], e.first));
					}
				}
			}
		}
	}
}

//set origin probabilities to default: uniform distribution
void traffic_network::set_origin_probabilities()
{
	origin_probabilities = std::vector<double>(number_of_nodes, 1.0 / number_of_nodes);
	random_origin = std::discrete_distribution<ULL>(origin_probabilities.begin(), origin_probabilities.end());
	recalc_mean_distances();
}

//set destination probabilities to default: uniform distribution
void traffic_network::set_destination_probabilities()
{
	destination_probabilities = std::vector<double>(number_of_nodes, 1.0 / number_of_nodes);
	random_origin = std::discrete_distribution<ULL>(destination_probabilities.begin(), destination_probabilities.end());
	recalc_mean_distances();
}

//set origin probabilities
void traffic_network::set_origin_probabilities(std::vector<double> param_origin_probabilities)
{
	assert(param_origin_probabilities.size() == number_of_nodes);

	origin_probabilities = param_origin_probabilities;

	//normalize probabilities
	double total_origin_prob = 0;
	for (ULL i = 0; i < number_of_nodes; ++i)
		total_origin_prob += origin_probabilities[i];
	for (ULL i = 0; i < number_of_nodes; ++i)
		origin_probabilities[i] /= total_origin_prob;

	random_origin = std::discrete_distribution<ULL>(origin_probabilities.begin(), origin_probabilities.end());
	recalc_mean_distances();
}

//set destination probabilities
void traffic_network::set_destination_probabilities(std::vector<double> param_dest_probabilities)
{
	assert(param_dest_probabilities.size() == number_of_nodes);

	destination_probabilities = param_dest_probabilities;

	//normalize probabilities
	double total_dest_prob = 0;
	for (ULL i = 0; i < number_of_nodes; ++i)
		total_dest_prob += destination_probabilities[i];
	for (ULL i = 0; i < number_of_nodes; ++i)
		destination_probabilities[i] /= total_dest_prob;

	random_destination = std::discrete_distribution<ULL>(destination_probabilities.begin(), destination_probabilities.end());
	recalc_mean_distances();
}

//compute mean distance with respect to the request distribution
void traffic_network::recalc_mean_distances()
{
	mean_pickup_distance = 0;
	mean_dropoff_distance = 0;

	// we need to exclude requests that want to go to the same node as where they already are.
	// this effectively redefines the normalization of probabilities for origin-destination pairs.
	// here we calculate these new probabilities.
	std::vector<std::vector<double>> redef_origin_dest_probs;
	for (ULL i = 0; i < number_of_nodes; ++i)
	{
		double total_redef_dest_prob = 0.0;
		redef_origin_dest_probs.push_back(std::vector<double>());
		for (ULL j = 0; j < number_of_nodes; ++j)
		{
			if (i != j)
			{
				total_redef_dest_prob += destination_probabilities[j];
			}
		}
		for (ULL j = 0; j < number_of_nodes; ++j)
		{
			if (i != j)
				redef_origin_dest_probs[i].push_back(origin_probabilities[i] * destination_probabilities[j] / total_redef_dest_prob);
			else
				redef_origin_dest_probs[i].push_back(0.0);
		}
	}
	for (ULL i = 0; i < number_of_nodes; ++i)
	{
		for (ULL j = 0; j < number_of_nodes; ++j)
		{
			mean_dropoff_distance += redef_origin_dest_probs[i][j] * get_network_distance(i, j);
			mean_pickup_distance += redef_origin_dest_probs[i][j] * get_network_distance(j, i);
		}
	}
}

//return distance in the network
double traffic_network::get_network_distance(ULL from, ULL to)
{
	return(network_distances[from][to]);
}

double traffic_network::get_p_2n()
{
	std::vector<int> p_l_to_go(number_of_nodes, 0);
	
	for (int i = 0; i < number_of_nodes; i++)
	{
		for (int j = 0; j < number_of_nodes; j++)
			p_l_to_go[(int)(get_network_distance(i, j) + 0.5)]++;
	}

	int sum1 = 0.0;
	int sum2 = 0.0;
	for (int i = 2; i < number_of_nodes; i++)
	{
		sum1 += p_l_to_go[i];
		for (int j = i; j < number_of_nodes; j++)
		{
			sum2 += p_l_to_go[j];
		}
	}

	//return (1.0 * p_l_to_go[1]) / (sum1 + p_l_to_go[1]);
	return (1.0 * sum1) / sum2;
}

//generate a new request based on the (uncorrelated) origin and destination probabilities
std::pair< ULL, ULL > traffic_network::generate_request()
{
	std::pair<ULL, ULL> request;

	request.first = random_origin(random_generator);
	request.second = request.first;

	while (request.second == request.first)
	{
		request.second = random_destination(random_generator);
	}

	return(request);
}

//return the shortest path from i to j in the form: r = ((i,t_i), (k_1,t_k_1), ... (k_n,t_k_n), (j,t_j))
//if i == j the route will have one node: r = ((i,t_i))
std::deque< std::pair<ULL, double> > traffic_network::find_shortest_path(ULL from, ULL to, double start_time, double velocity)
{
	route.clear();
	route.push_back(std::make_pair(from, start_time));

	double temp_route_time = -1;

	ULL current_node = from;
	double current_time = start_time;
	//pick nodes until the target is reached
	while (current_node != to)
	{
		//find next nearest node on the route
		temp_route_time = 1 + get_network_distance(current_node, to) / velocity;	//set distance to something larger than possible

		potential_route_nodes.clear();
		for (auto e : edgelist[current_node])
		{
			if (get_network_distance(e.first, to) / velocity + e.second / velocity < temp_route_time)
			{
				//if shorter distance found, clear the list of candidate nodes and add the node
				potential_route_nodes.clear();
				potential_route_nodes.push_back(e.first);
				temp_route_time = get_network_distance(e.first, to) / velocity + e.second / velocity;
			}
			else if (get_network_distance(e.first, to) / velocity + e.second / velocity == temp_route_time)
			{
				//add other node with the same distance to list of candidate nodes
				potential_route_nodes.push_back(e.first);
			}
		}
		//advance time along the route
		current_time += get_network_distance(current_node, potential_route_nodes[0]) / velocity;

																					   //choose randomly from all possible next nodes ( NOTE: this is NOT EXATCLY THE SAME(!) as choosing randomly from all possible routes, but good enough to randomize routes in regular topologies, e.g. a torus)
		current_node = potential_route_nodes[(ULL)(potential_route_nodes.size() * uniform01(random_generator))];

		//add the node to the route
		route.push_back(std::make_pair(current_node, current_time));
	}

	return(route);
}
