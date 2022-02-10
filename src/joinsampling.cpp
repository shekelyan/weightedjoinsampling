
typedef long double long_double;

//#include <external/robin_hood.h>

//#include <map>
//#include <set>


#define CSV_IO_NO_THREAD


namespace Global{
	const bool THREADS_ENABLED = false;
}

#include <array>

typedef std::array<uint8_t,16> offarray;

#include <joinsampling_db.h>

#include <approxdata/utils/utils.h>

class Input{

	shared_ptr<Database> db;
	
	map<string, string> parameters;
	
	string raw;
	
public:

	inline Input(int argc, char *argv[]){
	
		vector<string> args(argv + 1, argv + argc);
	
		stringstream ss;
	
		for (auto it = args.begin(); it != args.end(); it++)
			ss << " " << (*it);
		
		
		SQLParser p;
		p.parse(ss.str());
		
		for (auto it = p.expressions.begin(); it != p.expressions.end(); ++it){

			const Expression& e = (*it);		
			const ExpressionOperator& op = p.ops.at(e.op);
	
			if (UtilString::equals(op.name, "var")){
			
				const Expression* comment = p.getAncestor(e,"comment");
			
				if (comment){
		
					if (e.parent == comment->id){
					
						parameters[UtilString::toLowerCase(e.s)] = 1;
						
					} else {
		
						const Expression* eq = p.getAncestor(e,"equal");
				
						if (e.parent == eq->id)
							parameters[UtilString::toLowerCase(eq->childrenPtrs.at(0)->s)] = eq->childrenPtrs.at(1)->s;
					
					}
				}
			}
		}
		
		db = std::make_shared<Database>();
		
		
	
		Timer t1("load database "+parameters["db"]);
		
		t1.start();
		
		db->loadJson(parameters["db"], parameters);
		cout << endl;
		t1.end();
		cout << endl;
		
		Timer t2("execute query");
		t2.start();
		
		db->execute(ss.str(), parameters );
		
		cout << endl;
		t2.end();
		cout << endl;
		
	}

};


int main(int argc, char *argv[]){
	
	Input input(argc, argv);
	return 0;
	
}
