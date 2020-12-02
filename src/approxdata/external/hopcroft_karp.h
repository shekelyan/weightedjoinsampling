// C++ implementation of Hopcroft Karp algorithm for 
// maximum matching 

#include <list>
#include <queue>

using namespace std; 

// A class to represent Bipartite graph for Hopcroft 
// Karp implementation 
class BipartiteMatchingGraph 
{ 
	// m and n are number of vertices on left 
	// and right sides of Bipartite Graph 
	int m, n; 

	// adj[u] stores adjacents of left side 
	// vertex 'u'. The value of u ranges from 1 to m. 
	// 0 is used for dummy vertex 
	vector<vector<int>> adj; 

	// These are basically pointers to arrays needed 
	// for hopcroftKarp() 
	
	vector<int> pairU;
	vector<int> pairV;
	vector<int> dist;

public: 
	inline BipartiteMatchingGraph(int mm, int nn) : m(mm), n(nn), pairU(mm+1), pairV(nn+1), dist(mm+1){
		
		for (int j = 0; j <= mm; j++){
		
			vector<int> v;
			adj.push_back( v );
		}
	}
	void addEdge(int u, int v){
	
		adj[u].push_back(v); // Add u to v’s list. 
	}

	// Returns true if there is an augmenting path 
	inline bool bfs(){
	
		std::queue<int> Q; //an integer queue 

		// First layer of vertices (set distance as 0) 
		for (int u=1; u<=m; u++) 
		{ 
			// If this is a free vertex, add it to queue 
			if (pairU[u]==0) 
			{ 
				// u is not matched 
				dist[u] = 0; 
				Q.push(u); 
			} 

			// Else set distance as infinite so that this vertex 
			// is considered next time 
			else dist[u] = 99999999; 
		} 

		// Initialize distance to NIL as infinite 
		dist[0] = 99999999; 

		// Q is going to contain vertices of left side only. 
		while (!Q.empty()) 
		{ 
			// Dequeue a vertex 
			int u = Q.front(); 
			Q.pop(); 

			// If this node is not NIL and can provide a shorter path to NIL 
			if (dist[u] < dist[0]) 
			{ 
				// Get all adjacent vertices of the dequeued vertex u 
				
				for (auto it = adj[u].begin(); it !=adj[u].end(); it++) 
				{ 
					int v = *(it); 

					// If pair of v is not considered so far 
					// (v, pairV[V]) is not yet explored edge. 
					if (dist[pairV[v]] == 99999999) 
					{ 
						// Consider the pair and add it to queue 
						dist[pairV[v]] = dist[u] + 1; 
						Q.push(pairV[v]);
					} 
				} 
			} 
		} 

		// If we could come back to NIL using alternating path of distinct 
		// vertices then there is an augmenting path 
		return (dist[0] != 99999999); 

	}

	// Adds augmenting path if there is one beginning 
	// with u 
	bool dfs(int u){
	
		if (u == 0)
			return true;
	
		for (auto it = adj[u].begin(); it !=adj[u].end(); ++it){ 
			// Adjacent to u 
			int v = (*it); 

			// Follow the distances set by BFS 
			if (dist[pairV[v]] == dist[u]+1) 
			{ 
				// If dfs for pair of v also returns 
				// true 
				if (dfs(pairV[v]) == true) 
				{ 
					pairV[v] = u; 
					pairU[u] = v; 
					return true; 
				} 
			} 
		}

		// If there is no augmenting path beginning with u. 
		dist[u] = 99999999; 
		return false;
	}

	// Returns size of maximum matcing 
	inline int hopcroftKarp() 
	{

		// Initialize NIL as pair of all vertices 
		for (int u=0; u<=m; u++) 
			pairU[u] = 0; 
		for (int v=0; v<=n; v++) 
			pairV[v] = 0;
		
		// Initialize result 
		int result = 0; 

		// Keep updating the result while there is an 
		// augmenting path. 
		while (bfs()) 
		{ 
			// Find a free vertex 
			for (int u=1; u<=m; u++) 

				// If current vertex is free and there is 
				// an augmenting path from current vertex 
				if (pairU[u]==0){
				
					const int d = dfs(u);
			
			
					if (d){
					
						//cout << u << " " << pairU[u] << endl;
						result++; 
					}
					
				}
		} 
		return result; 
	}  
}; 


