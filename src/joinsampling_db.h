
#ifndef HEADER_JOINSAMPLING_DB
#define HEADER_JOINSAMPLING_DB

#include <joinsampling_lib.h>




class Database{

	vector<string> names;
	hash_map<string, shared_ptr<DatabaseTable> > tables;

	hash_map<string, string> views;
	
	bool debugcolnames = false;
	
	
public:	

	inline void loadJson(const string& file,const map<string, string>& parameters){
	
		nlohmann::json j;
		
		ifstream in(file);
		
		assert (in.is_open() );
	
		in >> j;
		
		in.close();
		
	
		assert (!j["tables"].empty() );
	
		nlohmann::json csv = j["tables"];
	
		Database* ret = this;
		
		for (json::iterator src = csv.begin(); src != csv.end(); ++src) {
	 
			const string key = src.key();
	 
			Tuple cols;
			Tuple sels;
			
			Tuple selnames;
			
			assert (!(*src)["csvcolumns"].empty() );
		
			nlohmann::json jcols = (*src)["csvcolumns"];
			
			nlohmann::json jsels = !(*src)["columns"].empty() ? (*src)["columns"] : jcols;
			
			char sepchar = ',';
			
			assert (!(*src)["csvfile"].empty() );
		
			const string file = (*src)["csvfile"];
		
			if (!(*src)["csvcomma"].empty()){
		
				const string s = (*src)["csvcomma"];
				sepchar = s[0];
			}
		
			bool csvheader = true;
		
			if (!(*src)["csvheader"].empty()){
		
				const bool s = (*src)["csvheader"];
				csvheader = s;
			}
			
			
			if (true){
			
				const string strcols = jcols;
				const string strsels = jsels;
				
				{
					const string cleaned = UtilString::cleanStringList(strcols);
				
					stringstream s_stream(cleaned);
					
					while(s_stream.good()) {
					  string substr;
					  getline(s_stream, substr, ',');
					  cols.push_back(substr);
					}
			    }
			    
			    {
					const string cleaned = UtilString::cleanStringList(strsels);
				
					stringstream s_stream(cleaned);
					
					while(s_stream.good()) {
					  string substr;
					  getline(s_stream, substr, ',');
					  sels.push_back(substr);
					}
			    }
			    
			} 
			
			if( (*src)["aliases"].empty() ){
			
				selnames = sels;
			
			} else {
			
				nlohmann::json jas = (*src)["aliases"];
		
				const string stras = jas;
				
				{
					const string cleaned = UtilString::cleanStringList(stras);
				
					stringstream s_stream(cleaned);
					
					while(s_stream.good()) {
					  string substr;
					  getline(s_stream, substr, ',');
					  selnames.push_back(substr);
					}
			    }
			
			}
			
			/*
			for (json::iterator col = jsels.begin(); col != jsels.end(); ++col){
			
				const string s = (*col);
				sels.push_back(s);
			}*/
			
			string f = file;
			
			f.reserve(100);
			
			for (auto it = parameters.begin(); it != parameters.end(); ++it){
			
				f = UtilString::replaced( f, "${"+it->first+"}", it->second);
			}
		
			auto tab = std::make_shared<CSVTable>(key, cols, sels, selnames, f, sepchar, false );
		
			if (!(*src)["csvdebug"].empty())
				tab->addColNames = (*src)["csvdebug"];
			
			if (debugcolnames)
				tab->addColNames = true;
		
			ret->addTable(key, tab);	
		}
		
		if (!j["views"].empty() ){
		
			nlohmann::json qs = j["views"];
			
			for (json::iterator qry = qs.begin(); qry != qs.end(); ++qry) {
				
				const string v = qry.value();
				
				addView(qry.key(), v);
			}
		
		}
		
	}

	
	inline Database& addTable(const string key2, shared_ptr<DatabaseTable> table){
		
		const string key = UtilString::toLowerCase(key2);
		
		tables[key] = table;
		names.push_back(key2);
		
		return (*this);
	}
	
	inline bool hasView(const string s) const{
		
		return views.count(UtilString::toLowerCase(s)) > 0;
	}
	
	inline bool hasTable(const string s) const{
		
		return tables.count(UtilString::toLowerCase(s)) > 0;
	}
	
	inline const string& getView(const string s) const{
		
		return views.at(UtilString::toLowerCase(s));
	}
	
	inline shared_ptr<DatabaseTable> getTable(const string& s) const{
		
		return tables.at(UtilString::toLowerCase(s));
	}
	
	inline Database& addView(const string key2, const string query){
		
		const string key = UtilString::toLowerCase(key2);
		
		if (views.count(key) > 0)
			views.at(key) = query;
		else
			views.insert( std::make_pair(key, query) );
		
		return (*this);
	}
	
	
	inline shared_ptr<JoinGraph> getJoinGraph(const string& cmd, const map<string, string>& parameters){
	
		std::set<string> processedTables;
		std::set<string> processedViews;
		Database* db = this;
	
		stringstream sqlstream;
	
		{
			SQLParser p;
		
			p.parse(cmd);
	
			for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){
	
				const Expression& e = (*it);		
				const ExpressionOperator& op = p.ops.at(e.op);
		
				if (UtilString::equals(op.name, "var")){
				
					const Expression* from = p.getAncestor(e,"from");
				
					if (from == NULL)
						from = p.getAncestor(e, "join");
				
					if (from){
			
						if (e.parent == from->id){
						
							const string cl = UtilString::cleanStringList(e.s, ' ', ' ');
							
							const string s = UtilString::getSplit(cl, " ", 0);
					
							if (db->hasView(s)){
							
								if (processedViews.count(s) == 0){
								
									sqlstream << "CREATE OR REPLACE VIEW " << s << " AS " << db->getView(s);
									
									processedViews.insert(s);	
								}	
							}
						}
					}
				}
			}
		}
		
		sqlstream << cmd;
		
		const string raw = sqlstream.str();
		
		auto query = std::make_shared<SQLExpression>(raw);
		
		{
		
			SQLParser& p = *query->parsed;
		
			for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){
	
				const Expression& e = (*it);		
				const ExpressionOperator& op = p.ops.at(e.op);
		
				if (UtilString::equals(op.name, "var")){
				
					const Expression* from = p.getAncestor(e, "from");
					
					if (from == NULL)
						from = p.getAncestor(e, "join");
					
					if (from){
			
						const string cl = UtilString::cleanStringList(e.s, ' ', ' ');
						const string s = UtilString::getSplit(cl, " ", 0);
						
						if (db->hasTable(s)){
						
							if (processedTables.count(s) == 0){
								
								query->addTable(s, db->getTable(s)->getColNames());
								processedTables.insert(s);
							}						
						}
					}
				}
			}
		}
		
		query->finalise();
		
		//query->print();
		
	
	
		
		SQLParser& p = *query->parsed;
		
		shared_ptr<JoinGraph> g = std::make_shared<JoinGraph>();
		
		
		g->setWeights(query);
	
		
		for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){

			const Expression& e = (*it);		
			const ExpressionOperator& op = p.ops.at(e.op);
	
			if (UtilString::equals(op.name, "var")){
			
				const Expression* from = p.getAncestor(e, "from");
				
				if (from == NULL)
					from = p.getAncestor(e, "join");
					
				if (from){
				
					const string cl = UtilString::cleanStringList(e.s, ' ', ' ');
					const string name = UtilString::getSplit(cl, " ", 0);
					
					const string alias = UtilString::getSplitLength(cl, " ") == 2 ? UtilString::getSplit(cl, " ", 1) : name;
					
					if (db->hasTable(name)){
					
						const Tuple& tabcols = db->getTable(name)->getColNames();
					
						Tuple newcols;
						Tuple selcols;
						Tuple joincols;
					
						//cout << "NAME " << name << " ALIAS " << alias << " COLS " << tabcols << endl;
					
						for (int j = 0; j < tabcols.size(); ++j){
						
							const string aliaskey = alias+"."+tabcols[j];
						
							if (query->selectedColumn(aliaskey)){
							
								selcols.push_back( tabcols[j] );
								newcols.push_back( query->renamedColumn(aliaskey) );
								
								
								if (query->joinColumn(aliaskey))
									joincols.push_back(query->renamedColumn(aliaskey));
							}
						
						}
						
						g->add( db->getTable(name)->projection( alias, selcols, newcols, false), joincols );
						
						//cout << name << " [" << tabcols << "] [" << selcols << "] " << alias << " ["<< newcols << "]" << "[" << joincols << "]" << endl;
						//Tuple selected;
						//Tuple joinCols;
						
						//Tuple 
					}
					
				}
			}
		}
		
		
		return g;
	}
	
	inline shared_ptr<WeightedRowSample> getSample(shared_ptr<JoinGraph> joingraph, JoinParameters& params, BasicAggregator& tableSizes, string seed = "abc"){
	
	
		
		shared_ptr<WeightedRowSample> sample;
		
		
	
		
		params.query = joingraph->parsed->parsed->repr();
		params.cyclic = joingraph->isCyclic();
		
		params.sampleSize = joingraph->parsed->limit;
		
		params.update();
		
		
		if (!params.weights){
			
			cout << endl << endl << "DISABLE WEIGHTS" << endl;
			joingraph->f->disable();
		}
		
		
		if (params.real("fk") > 0){
		
			params.maxSamples = params.sampleSize * params.real("fk");
		
			if (params.maxSamples > tableSizes.max()/2 )
				params.maxSamples = tableSizes.max()/2;
		
			params.print();
		
			BasicAggregator ag;
			
			auto& tables = joingraph->tables;
		
			for (long j = 0; j < tables.size(); ++j){
			
				ag.add(j, tables.at(j)->getSize()  );
			}
		
			const int root = ag.argMaxLast();
			
			auto in = tables.at(root);
			Tuple tin = in->getColNames();
			
			cout << "root table size " << in->getSize() << endl;
			
			auto smp = std::make_shared<WeightedRowSample>(tin, seed, params.maxSamples );	
	
			const bool weights = params.real("weights") != 0;
		
			{
				smp->initOnlineWOR();
				
				if (weights)
					in->startWeightScan(*joingraph->f, tin,tin);

				for (auto it = in->begin(); it != in->end(); ++it){
					smp->add(*it);
					
					in->weightRegister(*it);
				}
				
				if (weights)
					in->endWeightScan();
				
				smp->withReplacement();
			}
			
			Tuple tunion;
			
			if (weights)
			for (long jjj = 0; jjj < tables.size(); ++jjj){
			
				const long j = jjj == 0 ? root : jjj == root ? 0 : jjj;
			
				auto in1 = tables.at(j);
				Tuple tin1 = in1->getColNames();
				
				Tuple tin2 = UtilTuple::minus(tin1, tunion);
			
				if (tables.at(j)->getTotalWeight() == -1){
			
					in1->startWeightScan(*joingraph->f, tin1, tin2 );
				
					for (auto it = in1->begin(); it != in1->end(); ++it)
						in1->weightRegister(*it);
				
					in1->endWeightScan();
				}
				
				cout << tin2 << " max weight " << j << " " << in1->getMaxWeight() << endl;
				
				tunion = UtilTuple::merge(tunion, tin2);
			}
			
			
			Weight maxWeight = 1.0;
			
			for (long j = 0; j < tables.size(); ++j){
				
				assert (tables.at(j)->getMaxWeight()!= -1);
				
				maxWeight *= tables.at(j)->getMaxWeight();
			}
			
			cout << "max weight " << maxWeight << endl;
			
			
			
			long tabsSelected = 0;
			
			vector<shared_ptr<DatabaseTable>> tabs;
			
			
			Tuple tmerge;
			
			for (long k = 0; k < tables.size(); ++k ){
			
				for (long j = 0; j < tables.size(); ++j){
			
					if (tabsSelected & (1L << j) )
						continue;
						
					if (k == 0 && j != root)
						continue;
						
					if (tabs.size() == 0 || UtilTuple::intersection(tables.at(j)->getColNames(), tmerge).size() > 0){
					
						tabs.push_back( tables.at(j) );
						
						tmerge = UtilTuple::merge(tmerge, tables.at(j)->getColNames());
						tabsSelected |= 1L << j;
						break;
					}
				}
			}
			
			
			shared_ptr<DatabaseTable> sampletab = smp->getTable<DatabaseTable>(tables.at(root)->tableName);
			
			cout << "sampled root table size " << sampletab->getSize() << endl;
			
			
			for (long k = 1; k < tabs.size(); ++k ){
			
				
			
				DatabaseTable& a = *sampletab;
				DatabaseTable& b = *tabs.at(k);
				
				assert (a.inMemory() );
			
				Tuple ta = a.getColNames();
				Tuple tb = b.getColNames();
				
				cout << ta << " join " << tb << endl;
				
				Tuple tj = UtilTuple::intersection(ta, tb);
				
				assert (tj.size() > 0);
				
				Tuple tbnew = UtilTuple::minus(tb, ta);
								
				if (tbnew.size() == 0)
					continue;
				
				Tuple tc = UtilTuple::merge(ta, tbnew );
				
				vector<int> proj1 = UtilTuple::find(ta, tj);
				vector<long> hash1(proj1.size(), UtilHash::HASHMAX);
				
				vector<int> proj2 = UtilTuple::find(tb, tj);
				vector<long> hash2(proj2.size(), UtilHash::HASHMAX);
				
				vector<int> proj3 = UtilTuple::find(ta, ta);
				vector<int> proj4 = UtilTuple::find(tb, tbnew);
				
				RowHashMap<Row> hmap;
				
				for (auto it = a.begin(); it != a.end(); ++it){
				
					Row r = (*it);
				
					hmap.insert(r, proj1, hash1, r);
					
					
				}
					
				
				auto crows = std::make_shared<Rows>(tc.size() );	
				crows->reserve(a.getSize()*32);
				
				
				for (auto it = b.begin(); it != b.end(); ++it){
				
					const long longkey = hmap.getLongKey(*it, proj2, hash2);
				
					if (!hmap.hasKey(*it, proj2, hash2, longkey))
						continue;
				
					const Row& r1 = hmap.get(*it, proj2, hash2, longkey);
					const Row& r2 = *it;
					
					crows->push_back( r1, proj3, r2, proj4);
				}
				
				sampletab = std::make_shared<DatabaseTable>(a.tableName, tc, crows);
			}
			
			/*
			
			
			
			JoinGraph jg(*joingraph);
			
			cout << "root table size " << tables.at(root)->getSize() << endl;
			
			jg.tables.at(root) = smp->getTable<DatabaseTable>(tables.at(root)->tableName);
			
			cout << "sampled root table size " << jg.tables.at(root)->getSize() << endl;
			
			
			jg.shrinkCycles2(10000000000L, 0, false, subgraph(), subgraph(), true);
			
			
			
			assert (jg.tables.size() == 1);
			
			auto sampletab = jg.tables.at(0);
			
			cout << "joined table size " << sampletab->getSize() << endl;
			
			*/
			
			cout << "joined table size " << sampletab->getSize() << endl;
			
			Rows rows2( sampletab->getColNames().size() );
			
			rows2.reserve(sampletab->getSize() );
			
			Weighter wr(*joingraph->f, sampletab->getColNames());
			
			RandomNumbers rnd("test"+seed);
			
			long invalid = 0;
			
			for (auto it = sampletab->begin(); it != sampletab->end(); ++it){
				
				const Weight w = wr(*it);
				
				//cout << w << " " << maxWeight << endl;
				
				if (rnd.randomDouble() <= w/maxWeight){
					rows2.push_back(*it);
				} else {
				
					++invalid;
				}
			}
			
			sample = withoutrepl(rows2, sampletab->getColNames(), rows2.getSize(), seed);
			
			sample->setInvalidSize(invalid);
			sample->setMaxSize(params.maxSamples);
			sample->setPopulationWeight(tables.at(root)->getSize());
			sample->setPopulationSize(tables.at(root)->getSize());
			
			cout << sample->json() << endl;
			
			return sample;
		}
		
		
		Timer tapriori("phase0_apriori");
		
		tapriori.start();
		
		if (params.preload || params.buildIndex){
	
			if (params.preload || params.buildIndex){
	
				Timer t2("LOAD TABLES INTO MAIN MEMORY");
		
				t2.start();
		
				for (long j = 0; j < joingraph->tables.size(); ++j)
					joingraph->tables.at(j) = std::make_shared<DatabaseTable>(joingraph->tables.at(j));
				
				t2.end();
			}
		
			if (params.buildIndex){
		
				Timer t2("BUILD INDICES ON JOIN COLUMNS");
				t2.start();
		
				for (long j = 0; j < joingraph->tables.size(); ++j){
				
					Tuple tjoin;
					
					Tuple tt1 = joingraph->tables.at(j)->getColNames();
				
					for (long jj = 0; jj < joingraph->tables.size(); ++jj){
					
						if (jj == j)
							continue;
							
						Tuple tt2 = joingraph->tables.at(jj)->getColNames();
						
						tjoin = UtilTuple::merge(tjoin, UtilTuple::intersection(tt1, tt2));
					}
		
					cout << "tjoin " << tjoin << " for " << tt1 << endl;
				
					joingraph->tables.at(j)->buildIndices(tjoin);
					
					//cout << joingraph->tables.at(j)->toString(5) << endl;
				}
		
				t2.end();
			}
		}
		
		tapriori.end();
		
		
		// SIMPLIFY CYCLIC GRAPHS
		
		if (joingraph->isCyclic() && params.real("simplify") > 0 ){
		
			Timer tshrink("phase1_simplify");
			
			tshrink.start();
			
			subgraph black;
			subgraph white;
			
			double budget = params.real("simplify") * tableSizes.sum();
			
			joingraph->shrinkCycles(budget, false, black, white);
			
			
			
			cout << "join graph after joining along cycles: " << endl;
				
			cout << (joingraph->isCyclic() ?"CYCLIC " : "ACYCLIC ");
			cout << joingraph->toString() << endl;
			
			tshrink.end();
			
			
			if (joingraph->isCyclic() )
				params.buildIndex = true;
			else
				params.cyclic = false;
				
			params.update();
				
			
			if (params.buildIndex){
	
				Timer t2("BUILD INDICES");
				t2.start();
			
				for (long j = 0; j < joingraph->tables.size(); ++j)
					joingraph->tables.at(j)->buildIndices();
			
				t2.end();
			}
			
		}
		cout << "A1" << endl;
		
		
		if (params.real("inversion") > 0){
		
			Timer t("inversion");
			t.start();
			
			Timer tshrink("phase1_simplify");
			
			tshrink.start();
		
			joingraph->shrinkCycles2(10000000000L, 1, false, subgraph(), subgraph(), true);
			
			tshrink.end();
			
			WeightingFunction& f = *joingraph->f;
			
			cout << joingraph->toString() << endl;
			
			DatabaseTables tabs(joingraph->tables);
			
			//cout << endl;
			//cout << "JOIN SIZE " << tabs.joinsize() << endl;
			//cout << endl;
			
			cout << "A2" << endl;
			
			
			Timer tdraw("phase2_sampling");
			
			tdraw.start();
			
			cout << "A3" << endl;
			
			//auto smp = tabs.sample(params.sampleSize, f, "abc123");
			
			sample = tabs.sample(params.sampleSize, f, "abc123");
			
			tdraw.end();
			
			//const Tuple& smpCols = smp->getColNames();
			
			/*
			sample = std::make_shared<WeightedRowSample>(smpCols, params.get("seed", "abc123"), smp->getSize() );
			
			
			
			sample->setPopulationWeight();
			sample->setPopulationSize();*/
			
			t.end();
			
			cout << sample->json() << endl;
			
			return sample;
		}
		
		
		// ACYCLIC GRAPHS
		
		if (true){
		
			Timer tdraw("phase2_sample");
			
			tdraw.start();
		
			cout << "A4" << endl;
			
			params.accrate = 1.0;
			
			subgraph cyclebreaks = joingraph->getCyclebreaks(params.accrate, params.real("break", 0));
			
			subgraph semi;
			
			params.update();
			
			const bool breakcycles = cyclebreaks.size() > 0;
			
			if ((params.real("semijoin") >= 1) &&  breakcycles){
			
				cout << joingraph->toString() << endl;
				cout << "without semijoins " << joingraph->toString(cyclebreaks) << endl;
			
				for (auto it = cyclebreaks.begin(); it != cyclebreaks.end(); ++it){
				
					const int e = *it;
					
					if (joingraph->g.isEdge(e) ){
					
						auto ee = joingraph->g.getConst(e);
						
						const int c = joingraph->g.pos(e);
						
						Tuple t = UtilTuple::tuple( joingraph->edgeLabels.at(c) );
						
						const int a = joingraph->g.pos(ee.at(2));
						const int b = joingraph->g.pos(ee.at(1));
						 
						 
						const int e1 = joingraph->addLink(a, joingraph->tables.at(b), t);
						
						semi.insert(e1);
						
						if (params.real("semijoin") >= 3){
						
							vector<int> sp = joingraph->g.shortestPath(b,a, true, true, cyclebreaks);
							
							for (int k = 0; k < (sp.size()-1)/2; ++k){
						
								int prev = joingraph->tables.size()-1;
							
								int e3 = sp[k*2+1];
								int v3 = sp[k*2+2];
							
								Tuple t3 = UtilTuple::tuple( joingraph->edgeLabels.at( joingraph->g.pos(e3) ) );
							
								const int e4 = joingraph->addLink(prev, joingraph->tables.at(joingraph->g.pos(v3) ), t3);
							
								semi.insert(e4);
							}
						}
						
						if (params.real("semijoin") >= 2){
						
							const int e2 = joingraph->addLink(b, joingraph->tables.at(a), t);
							semi.insert(e2);
						
							if (params.real("semijoin") >= 4){
						
								vector<int> sp = joingraph->g.shortestPath(a,b, true, true, cyclebreaks);
						
								for (int k = 0; k < (sp.size()-1)/2; ++k){
						
									int prev = joingraph->tables.size()-1;
							
									int e3 = sp[k*2+1];
									int v3 = sp[k*2+2];
								
									Tuple t3 = UtilTuple::tuple( joingraph->edgeLabels.at( joingraph->g.pos(e3) ) );
							
									const int e4 = joingraph->addLink(prev, joingraph->tables.at(joingraph->g.pos(v3) ), t3);
							
									semi.insert(e4);
								}
							}
						}
						
						
						
						//cout << "[" << joingraph->label(a) << "]--" << joingraph->label(e) << "--["<< joingraph->label(b)<< "]" << endl;
						
						
					}
				}
				
				cout << "with semijoins " << joingraph->toString(cyclebreaks, semi) << endl;
			}
			
			cout << "A5" << endl;
			
			subgraph fk = params.real("fk") == 2 ? joingraph->g.component(joingraph->g.vertexIdFromPos(0), false, true) : subgraph();
			
			cout << "foreign key edges " << fk.size() << endl;
			
			auto joinTree = joingraph->getJoinTree(breakcycles, cyclebreaks, semi, fk);
			
			
			cout << joinTree->toString() << endl;
			
			cout << "A6" << endl;
			
			//joingraph->f->print();
			
			Timer t4("SAMPLE");
			t4.start();
		
			sample = joinTree->getSample(params, params.get("seed", "abc123"));
		
			cout << endl;
			t4.end();
			
			tdraw.end();
			
			cout << "A7" << endl;
			
			cout << sample->json() << endl;
			
			return sample;
			
		}
	
	
	
	
	
	}
	
	inline void execute(const string& cmd, const map<string, string>& parameters){
	
		
		Timer t1("JOIN GRAPH");
		t1.start();
		
		auto joingraph = getJoinGraph(cmd, parameters);
		
		
		BasicAggregator tableSizes;
		 
		{
			Timer t2("phase0_sizes");
			t2.start();
	
			for (long j = 0; j < joingraph->tables.size(); ++j){
				
				joingraph->tables.at(j)->init();
				
				tableSizes.add(joingraph->tables.at(j)->getSize() );
			}
			t2.end();
		}
		
		JoinParameters params(parameters);
		
		cout << endl << "SEED " << params.get("seed") << endl << endl;
		
		
		const double sc = params.real("smp", 1.0);
		
		cout << "inputted join graph: " << endl;
		cout << joingraph->toString() << endl;
		
		
		for (long j = 0; j < joingraph->tables.size(); ++j){
		
			joingraph->tables.at(j)->enableSampling(sc*tableSizes.max() );
			
			cout << joingraph->tables.at(j)->tableName << " " << joingraph->tables.at(j)->getColNames() << " " << joingraph->tables.at(j)->getSize() << endl;
		}
		
		t1.end();
		
		
		shared_ptr<WeightedRowSample> sample;
		
		for (int jj = 0; jj < 100; ++jj){
		
			auto samp1 = getSample(joingraph, params, tableSizes, "abc"+to_string(jj) );
		
			if (sample){
				
				sample->push_back(*samp1);
		
			} else {
		
				sample = samp1;
			}
			
			if (sample->getSize() >= params.sampleSize)
				break;
		}
		
		if (sample){
		
		} else {
		
			assert (false);
		}
		
		
		cout << endl;
		cout << endl;
		cout << endl;
		
		cout << sample->json() << endl;
		
		
		joingraph->f->verbose = false;
		
		Weighter w(*(joingraph->f), sample->cols);
		
		long k = 10;
		
		cout << endl << "First ten samples: " << endl;
		
		for (auto it = sample->begin(); it != sample->end(); ++it){
		
			cout << (*it) << " weight " << w(*it) << endl;
			
			--k;
			if (k == 0)
				break;
		}
		
		cout << endl << "Last ten samples: " << endl;
		
		k = sample->getSize();
		
		for (auto it = sample->begin(); it != sample->end(); ++it){
		
			if (k < 10)
				cout << (*it) << " weight " << w(*it) << endl;
			
			--k;
		}
		
		if (params.defined("write")){
		
			stringstream ss;
			
			for (int j = 0; j < sample->cols.size(); ++j)
				ss << sample->cols[j] << "_";
			
			ss << params.get("write");
		
			const string filename = ss.str();
		
			ofstream file (filename);
			
			if (file.is_open()){
			
				const char sep = UtilString::endsWith( UtilString::toLowerCase(filename), ".tbl") ? '|' : ',';
				
				for (auto it = sample->begin(); it != sample->end(); ++it){
				
					const Row& r = (*it);
					
					file << r.get(0);
					
					for (int i = 1; i < r.size(); ++i)
						file << sep << r.get(i);
					
					file << endl;
				}
				
				file.close();
			  }
		}
		
		if (params.defined("ks") && params.enabled("ks") ){
		
			
			if (sample){
				/**/
				
				for (int k = 0; k < 2; k++){
					
					if (k == 0 && sc == 1.0)
						continue;
					
					auto joingraph = getJoinGraph(cmd, parameters);
	
					BasicAggregator tableSizes2;
	
					{
						Timer t2("GET TABLE SIZES");
						t2.start();

						for (long j = 0; j < joingraph->tables.size(); ++j){
			
							joingraph->tables.at(j)->init();
							
							cout << joingraph->tables.at(j)->tableName << " " << joingraph->tables.at(j)->getColNames() << " " << joingraph->tables.at(j)->getSize() << endl;
			
							tableSizes2.add(joingraph->tables.at(j)->getSize() );
						}
						t2.end();
					}
					
					double smp = k == 1 ? 1.0 : sc;
					
					if (smp < 1.0){
					
						for (long j = 0; j < joingraph->tables.size(); ++j){
						
							joingraph->tables.at(j)->enableSampling(smp*tableSizes2.max() );
						}
					}
					
					Timer t("validate join-sample");
					t.start();
				
					joingraph->shrinkCycles2(10000000000L, 1, false, subgraph(), subgraph(), true);
				
					WeightingFunction& f = *joingraph->f;
				
					cout << "tabs" << endl;
				
					DatabaseTables tabs(joingraph->tables);
					
					/*
					const long double jw = tabs.joinweight(f);
			
					cout << "smp " << smp << " ks.joinweight " << jw << endl; 
					
					const long js = tabs.joinsize();
					
					cout << "smp " << smp << " ks.joinsize " << js << endl;
					*/
					
					auto dbt = sample->getTable<DatabaseTable>("sample");
					
					cout << "DBT1 " << dbt->getSize() << endl;
					
					Tuple joincols = tabs.sj().joinCols;
					
					for (int i = 0; i < joincols.size(); ++i){
					
						if (UtilTuple::find(dbt->getColNames(), joincols[i]) == -1){
						
							if (UtilString::contains(joincols[i], "_")){
							
								dbt = dbt->addColumnUnion(UtilTuple::tuple(joincols[i] ));
								cout << "added " << joincols[i] << endl;
							} else {
							
								assert (false);
							}
						}
					}
					
					dbt->enableSampling(params.sampleSize);
					
					dbt = std::make_shared<DatabaseTable>(dbt);
					
					cout << "DBT2 " << dbt->getSize() << endl;
					
					for(int i = 0; i < 1; ++i){
					
						cout << "dstat " << endl;
					
						const double dst = tabs.dstat(*dbt, f);
						
						
						if ((k == 0) || (sc == 1.0) ){
						
							cout << "ks_smp_pop.dstat" << " " << dst << endl;
						}
						
						if ((k == 1) && (sc < smp)){
						
							cout << "ks_smp_smp.dstat" << " " << dst << endl;
						}
						
						if (sc < smp){
							
							cout << "ks_approx.dstat" << " " << dst << endl;
							
						} else {
						
							cout << "ks_exact.dstat" << " " << dst << endl;	
						}
						
					}
					
					t.end();
				
				}
				
				return;
			}
		}
	
		
		//cout << sample->toString(10) << endl;
	}
};



#endif
