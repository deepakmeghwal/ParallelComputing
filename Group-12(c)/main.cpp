/**
	@ Abhinav Bollam  - 150102075
	@ Saurabh Dhall   - 150102060
	@ Deepak Meghwal  - 150108009

	Implementation of Distributed Parallel Cooperative Coevolution based Multi-objective
	optimization algorithm for large scale optimization problems:

 	https://ieeexplore.ieee.org/abstract/document/7867084
 */
#include <bits/stdc++.h>
#include <time.h>
using namespace std;

int num_Of_Variables;
int num_Of_Objective_Functions;
int num_Of_Samples, populationSize;
int max_num_Of_Evaluations =  50000;
int evals = 0;
int mode  = 1;
double e;
double pi = 3.141592;
double t =20;
double delta = 0.9;
double cr = 1.0;
double F =  0.5;
double mutp;
double muteta = 20.0;
double errorlimit = 0.0001;
int num_changes = 2;




////////////////////////////////////////////////////////////////////////////////////////
//*********************************** DTLZ Problem sets ******************************//
////////////////////////////////////////////////////////////////////////////////////////


vector<double> evaluate(vector<double> & var){
	int k = num_Of_Variables - num_Of_Objective_Functions +1;

	vector<double> function_Values(num_Of_Objective_Functions);
	vector<double> angles(num_Of_Objective_Functions - 1);
	
	double h = 0.0;
	double exp = 100.0;

	switch(mode) {

		/////////////////////////// DTLZ1 /////////////////////////

		case 1:
			for (int i = num_Of_Variables - k; i<num_Of_Variables; i++)
				h += ( (var[i] - 0.5)*(var[i] - 0.5) - cos(20.0 * pi * (var[i] - 0.5)) );

			h = 100 * ( k + h );
			for (int i = 0; i < num_Of_Objective_Functions; i++)
				function_Values[i] = (1.0 + h) * 0.5 ;

			for (int i = 0; i < num_Of_Objective_Functions; i++){
				for (int j = 0; j < num_Of_Objective_Functions - (i + 1); j++)
					function_Values[i] *= var[j];
				if (i != 0){
					int temp = num_Of_Objective_Functions - (i + 1);
					function_Values[i] *= 1 - var[temp];
				} 
			}
			break;
		
		////////////////////////// DTLZ2 //////////////////////////

		case 2: 
			for (int i = num_Of_Variables - k; i < num_Of_Variables; i++)
    			h += (var[i] - 0.5)*(var[i] - 0.5);

  			for (int i = 0; i < num_Of_Objective_Functions; i++)
    			function_Values[i] = 1.0 + h;

			for (int i = 0; i < num_Of_Objective_Functions; i++){
				for (int j = 0; j < num_Of_Objective_Functions - (i + 1); j++)
					function_Values[i] *= cos(var[j]*0.5*pi);
				if (i != 0){
					int temp = num_Of_Objective_Functions - (i + 1);
					function_Values[i] *= sin(var[temp]*0.5*pi);
				} 
			}
			break;

		//////////////////////// DTLZ3 ///////////////////////////

		case 3:
			for (int i = num_Of_Variables - k; i<num_Of_Variables; i++)
				h += ( (var[i] - 0.5)*(var[i] - 0.5) - cos(20.0 * pi * (var[i] - 0.5)) );

			h = 100.0 * (k + h);
			for (int i = 0; i < num_Of_Objective_Functions; i++)
    			function_Values[i] = 1.0 + h;

    		for (int i = 0; i < num_Of_Objective_Functions; i++){
				for (int j = 0; j < num_Of_Objective_Functions - (i + 1); j++)
					function_Values[i] *= cos(var[j]*0.5*pi);
				if (i != 0){
					int temp = num_Of_Objective_Functions - (i + 1);
					function_Values[i] *= sin(var[temp]*0.5*pi);
				} 
			}
			break;

		////////////////////// DTLZ4 //////////////////////////////

		case 4:
			for (int i = num_Of_Variables - k; i < num_Of_Variables; i++)
    			h += (var[i] - 0.5)*(var[i] - 0.5);

    		for (int i = 0; i < num_Of_Objective_Functions; i++)
    			function_Values[i] = 1.0 + h;

    		for (int i = 0; i < num_Of_Objective_Functions; i++){
				for (int j = 0; j < num_Of_Objective_Functions - (i + 1); j++)
					function_Values[i] *= cos( pow(var[j],exp) * 0.5 * pi);
				if (i != 0){
					int temp = num_Of_Objective_Functions - (i + 1);
					function_Values[i] *= sin( pow(var[temp],exp) * 0.5 * pi);
				} 
			}
			break;

		////////////////////// DTLZ5 ///////////////////////////////

		case 5:
			for (int i = num_Of_Variables - k; i < num_Of_Variables; i++)
    			h += (var[i] - 0.5)*(var[i] - 0.5);

    		e = pi / (4.0 * (1.0 + h));

    		angles[0] = var[0] * 0.5 * pi;
    		for (int i=0; i < (num_Of_Objective_Functions -1); i++)
    			angles[i] = e * (1.0 + 2.0 * h * var[i]);

    		for (int i = 0; i < num_Of_Objective_Functions; i++)
    			function_Values[i] = 1.0 + h;

    		for (int i = 0; i < num_Of_Objective_Functions; i++){
				for (int j = 0; j < num_Of_Objective_Functions - (i + 1); j++)
					function_Values[i] *= cos( angles[j] );
				if (i != 0){
					int temp = num_Of_Objective_Functions - (i + 1);
					function_Values[i] *= sin( angles[temp] );
				} 
			}
			break;

		///////////////////// DTLZ6 ////////////////////////////////

		case 6:
			for (int i = num_Of_Variables - k; i < num_Of_Variables; i++)
    			h += pow(var[i], 0.1);

    		e = pi / (4.0 * (1.0 + h));			

    		angles[0] = var[0] * 0.5 * pi;
    		for (int i=0; i < (num_Of_Objective_Functions -1); i++)
    			angles[i] = e * (1.0 + 2.0 * h * var[i]);

    		for (int i = 0; i < num_Of_Objective_Functions; i++)
    			function_Values[i] = 1.0 + h;

    		for (int i = 0; i < num_Of_Objective_Functions; i++){
				for (int j = 0; j < num_Of_Objective_Functions - (i + 1); j++)
					function_Values[i] *= cos( angles[j] );
				if (i != 0){
					int temp = num_Of_Objective_Functions - (i + 1);
					function_Values[i] *= sin( angles[temp] );
				} 
			}
			break;

		////////////////////////// DTLZ7 /////////////////////////////

		default :
			for (int i = num_Of_Variables - k; i < num_Of_Variables; i++)
    			h += var[i];

    		h = 1 + (9.0 * h)/k;

    		for (int i=0; i < (num_Of_Objective_Functions -1); i++)
    			function_Values[i] = var[i];

    		e = 0.0;
    		for (int i=0; i < (num_Of_Objective_Functions -1); i++)
    			e += ( function_Values[i]/(1.0+h) ) * (1 + sin(3.0 * pi * function_Values[i]) );

    		e = num_Of_Objective_Functions - e;
    		function_Values[num_Of_Objective_Functions - 1] = (1+h)*e;

    		break;
    	}

    

	return function_Values;

}




vector<double> createvariables(){
	vector<double> variable_ (num_Of_Variables);
	for ( int i=0; i<num_Of_Variables ;  i++){
		variable_[i] = rand()/(double)(RAND_MAX);
	}
	return variable_;
}


////////////////////////////////////////////////////////////////////////////////////////
//*********************************** utility functions ********************************
////////////////////////////////////////////////////////////////////////////////////////


double calError(vector< vector<double> > & A, vector< vector<double> > & B){
	double sum = 0;
	for (int w=0; w<num_Of_Samples; w++){
		for (int r=0; r<num_Of_Objective_Functions; r++){
			sum += pow(abs(A[w][r]-B[w][r]) , 2);
		}
	}
	return sqrt(sum);
}

double distance(vector<double> a, vector<double> b){
	 double ans =0;
	 for(int i=0;i<a.size();i++){
	 	ans += (a[i]-b[i])*(a[i]-b[i]); 
	 }
	 return sqrt(ans);
}

void smallsort(vector<double> & val, vector<double> & ind){
	for (int i=0;i<t;i++){
		 for(int j=i+1; j<num_Of_Samples; j++){
		 	if(val[i]>val[j]){
		 		double temp= val[i]; val[i]= val[j]; val[j] = temp;
		 			   temp = ind[i]; ind[i]= ind[j]; ind[j] = temp;
		 	}
		 }
	}
}


void shuffle(vector<double> & perm){
	for(int i= perm.size()-1;i>0;i--){
		int j= rand()%(i+1);
		double temp = perm[i]; perm[i]= perm[j]; perm[j] = temp;
	}
}

int randint(int l, int u){
	return (rand()%(u-l+1))+l;
}

/////////////  approximate difference between a, z  relative to b ////////////////
double approximate(vector<double> a, vector<double> b, vector<double> z){
	double optimal_diff = 0.0;
	double maxval = -1*pow(10,30);
	
	for( int w=0; w<num_Of_Objective_Functions; w++){
		double diff = fabs( a[w] - z[w]);
		double rel_diff ;
		if(b[w] == 0){
			rel_diff = 0.0001 * b[w];
		}else{
			rel_diff = diff * b[w];
		}
		if(rel_diff > maxval){
			maxval = rel_diff;
		}
	}
	optimal_diff = maxval;
	return optimal_diff;
}




////////////////////////////////////////////////////////////////////////////////////////
//*********************************** main ***********************************************
////////////////////////////////////////////////////////////////////////////////////////

int main(int arg, char** argv ){
	
	cout<<"DPCCMOEA - LSOP is implemented for Multiobjective problem set for Real Values."<<endl;
	cout<<"Choose type of problem set.\n 	1. DTLZ1 \n 	2. DTLZ2 \n 	3. DTLZ3 \n 	4. DTLZ4 \n"<<
	 	" 	5. DTLZ5 \n 	6. DTLZ6 \n 	7. DTLZ7  "  << endl;

	cout<<"Enter a number [1-7]: "<<endl;
	cin>> mode;

	cout<<"Enter number of Variables :  "<<endl;
	cin>> num_Of_Variables;

	cout<<"Enter number of Objective Functions from set[ 2, 3, 5 ]: "<<endl;
	cin>> num_Of_Objective_Functions;

	cout<<"The maximum size of population is 1000. Enter size of population : "<<endl;
	cin>> num_Of_Samples;
	populationSize = num_Of_Samples;

	cout<<num_Of_Samples<<" samples with each "<<num_Of_Variables<<" variables are randomly genearated for DPCCMOEA-LSOP "<<endl;

	
	////////////////////////////////////////////////////////////////////////////////////
	//////////////////// MOEAD implementation //////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////



	////////////  Reading data///////

	vector<vector<double> > lambda(num_Of_Samples);
	for(int i=0;i<num_Of_Samples; i++){
		lambda[i] = vector<double>(num_Of_Objective_Functions);
	}

	ostringstream os;
    os << "data/Weight/" << "W" << num_Of_Objective_Functions << "D_" << to_string(1000) << ".dat";
    string FileName;
    FileName = os.str();
    //cout<<FileName<<endl;

    ifstream in(FileName.c_str());
	if( !in ) {
      cout << "Failed when reading from file: : " <<FileName << endl;
      exit(-1);
    }
    else{
	    /// Reading file
	    int i = 0;
	    int j = 0;
	    string str;
	    while (getline(in, str)) {
	    	istringstream iss(str);
	    	j = 0;
	      
	    	while (iss) {
	        	string token;
	        	iss >> token;
	        	if (token.compare("")!=0) {
	          		double value = atof(token.c_str());
	         		lambda[i][j] = value;
	          		j++;
	        	} 
	      	} // while
	      	i++;
	      	if(i==num_Of_Samples)
	      		break;
	    } // while
	    in.close();
	} // else

	cout<<"Reading data successful"<<endl;

	/////////////  Start of DPCCMOEA Algo /////////////////////
	int start = clock();

	//////////  Finding nearby variables using DVA //////////// 

	vector<double> Neigh(num_Of_Samples);
	vector<double> Neighind(num_Of_Samples);

	vector<vector<double> > ansNeigh(num_Of_Samples);
	for(int i=0;i<num_Of_Samples; i++){
		ansNeigh[i] = vector<double>(t);
	}

	for(int i=0; i<num_Of_Samples; i++){
		for(int j=0; j<num_Of_Samples; j++){
			Neigh[j] = distance(lambda[i], lambda[j]);
			Neighind[j] = j; 
		}
		smallsort(Neigh, Neighind);

		for(int j=0;j<t;j++){
			ansNeigh[i][j] = Neighind[j];
		}
	}

	cout<<" Finding Nearby variables successful"<<endl;

	////////////   Initialize Populations //////////

	vector<vector<double> > ansPopVar(num_Of_Samples);
	for(int i=0;i<num_Of_Samples; i++){
		ansPopVar[i] = vector<double>(num_Of_Variables);
	}

	vector<vector<double> > ansPopFun(num_Of_Samples);
	for(int i=0;i<num_Of_Samples; i++){
		ansPopFun[i] = vector<double>(num_Of_Objective_Functions);
	}

	vector<vector<double> > prevAnsPopFun(num_Of_Samples);
	for(int i=0;i<num_Of_Samples; i++){
		prevAnsPopFun[i] = vector<double>(num_Of_Objective_Functions);
	}

	for(int i=0;i<num_Of_Samples; i++){
		ansPopVar[i] = createvariables();
		ansPopFun[i] = evaluate(ansPopVar[i]);
		evals++;
	}

	cout<<" Population initialization successful"<<endl;

	/////////// Find overall minimum function values //////

	vector<double> funMinVal(num_Of_Objective_Functions);
	for(int i=0;i<num_Of_Objective_Functions;i++){
		funMinVal[i] = pow(10,30);
	}

	for(int i=0;i<num_Of_Samples;i++){
		for(int j=0;j<num_Of_Objective_Functions;j++){
			 if(funMinVal[j] > ansPopFun[i][j]){
			 	funMinVal[j] = ansPopFun[i][j];
			 }
		}
	}

	cout<<" Finding overall minimum initiaited"<<endl;

	//////////  Population evolution  ///////////////////

	cout<<" Started population evolution"<<endl;

	while(evals< max_num_Of_Evaluations){
		vector<double> perm(num_Of_Samples);
		for(int i=0;i<num_Of_Samples;i++){
			perm[i]= i;
		}
		shuffle(perm);
		for(int i=0; i<num_Of_Samples; i++){
			int curr_id = perm[i];
			int type;
			double rval = rand()/(double)(RAND_MAX);

			//  choosing type : Neighbourhood or whole set
			if( rval < delta){
				type = 1;
			}else{
				type = 2;
			}

			// Optimization of the population using Differential Evolution Crossover
			vector<int> p;						// select 2 parents for DE of current Sample. 
			while( p.size() < 2){
				int temp;
				if(type == 1){
					int r = randint(0, t-1);
					temp = ansNeigh[curr_id][r];
				}else{
					temp = randint(0, num_Of_Samples-1);
				}

				bool flag = 1;
				for( int w =0; w<p.size(); w++){
					if(p[w]==temp){
						flag=0;
						break;
					}
				}
				if(flag){
					p.push_back(temp);
				}
			}

			// DE crossover
			vector<double> childvar(num_Of_Variables);
			double rval2 =  rand()/(double)(RAND_MAX);
			double rval3 =  rand()/(double)(RAND_MAX);
			for( int w=0; w<num_Of_Variables ; w++){
				childvar[w] = ansPopVar[curr_id][w] + F*( ansPopVar[p[0]][w] - ansPopVar[p[1]][w]);
				if(childvar[w] < 0.0){
					childvar[w] = 0.0;
				}else if(childvar[w] > 1.0){
					childvar[w] = 1.0;
				}
			}
			if(type == 2){
				if(rval2 <= 0.5){
					for( int w=0; w<num_Of_Variables ; w++){
						childvar[w] = ansPopVar[curr_id][w];
					}
				}else{
					if(rval3 <= 0.5){
						for( int w=0; w<num_Of_Variables ; w++){
							childvar[w] = ansPopVar[p[0]][w];
						}
					}else{
						for( int w=0; w<num_Of_Variables ; w++){
							childvar[w] = ansPopVar[p[1]][w];
						}
					}
				}
			}

			// Polynomial mutation
			for(int w=0; w<num_Of_Variables; w++){
				double rval4 = rand()/(double)(RAND_MAX);
				mutp = 1/num_Of_Variables;
				if(rval4 <= mutp){
					double xmin = 0.0, xmax= 1.0;
					double delta1 = (childvar[w] - xmin)/(xmax - xmin);
					double delta2 = (xmax - childvar[w])/(xmax - xmin);
					double rval5  = rand()/(double)(RAND_MAX);
					double coeff, deltaq, mid, base;
					if(rval5 <= 0.5){
						coeff = (1.0 - 2*rval5);
						mid   = pow(1-delta1, muteta+1 );
						base  = (2*rval5 + base); 
						deltaq= pow(base, muteta+1 ) - 1.0; 

					}else{
						coeff = (2*rval5 - 1.0);
						mid   = pow(1-delta2, muteta+1 );
						base  = (2*(1-rval5) + base );
						deltaq= 1.0 - pow(base, muteta+1);
					} 

					childvar[w] = childvar[w] + deltaq*(xmax-xmin);
					if(childvar[w] < xmin){
						childvar[w] = xmin;
					}else if( childvar[w] > xmax){
						childvar[w] = xmax;
					}
				}
			}

			////////// Evaluate Child vector and update overall min/////////
			vector<double> childfun = evaluate(childvar);

			for(int w=0; w<num_Of_Objective_Functions; w++){
				 if(funMinVal[w] > childfun[w]){
				 	funMinVal[w] = childfun[w];
				 }
			}

			////////// Update Solutions  ////////////////////////////////

			int changes = 0;
			int siz = 0;
			type==1 ? siz=(int)t : siz=(int)num_Of_Samples;
			
			vector<double> rot(siz);
			for(int w=0; w<siz; w++){
				rot[w] =  w;
			}
			shuffle(rot);

			for(int w=0; w<siz; w++){
	
				int q;
				if( type==1 ){
					q = ansNeigh[curr_id][rot[w]];
				}else{
					q = rot[w];
				}
				
				double app1, app2;
				app1 = approximate(ansPopFun[q], lambda[q], funMinVal);
				app2 = approximate(childfun    , lambda[q], funMinVal);

				if(app2 < app1){
					ansPopVar[q] = childvar;
					ansPopFun[q] = evaluate(childvar);
					
					changes++; 
				}

				if(changes >= num_changes){
					break;
				}

			}


		}

		


		if( evals % 200 == 0){
			for (int w=0; w<num_Of_Samples; w++){
				for (int r=0; r<num_Of_Objective_Functions; r++){
					prevAnsPopFun[w][r] = ansPopFun[w][r];
				}
			}
		}

		if( evals % 200 == 1 && evals > 1 ) {
			double error = calError(ansPopFun, prevAnsPopFun);

			if(error < errorlimit){
				break;
			}
		}

		evals++;

		if(evals%200 == 0)
			cout<< "Number of evaluations so far "<<evals<<endl;
		
	}


	
	int end   = clock();
  	double cycles = (double) (end - start);
  	double secs = cycles / CLOCKS_PER_SEC;


	ofstream op_file;
	op_file.open("Var_DTLZ"+to_string(mode)+"_"+to_string(num_Of_Objective_Functions)+".txt");

	for(int w=0; w<num_Of_Samples; w++){
		for(int q=0;q< num_Of_Variables; q++){
			op_file<<ansPopVar[w][q];
			if( q < (num_Of_Variables-1) ){
				op_file<<" ";
			}
		}
		op_file<<"\n";
	}

	op_file.close();

	ofstream opfun_file;
	opfun_file.open("Fun_DTLZ"+to_string(mode)+"_"+to_string(num_Of_Objective_Functions)+".txt");

	for(int w=0; w<num_Of_Samples; w++){
		vector<double> objfun = evaluate(ansPopVar[w]);
		for(int q=0;q< num_Of_Objective_Functions; q++){
			opfun_file<<objfun[q];
			if( q < (num_Of_Objective_Functions-1) ){
				opfun_file<<" ";
			}
		}
		opfun_file<<"\n";
	}

	opfun_file.close();

	////////////////  Results ////////////////////////

	 // Result messages
	cout << endl;
	cout << "Type of Problem set choosen - DTLZ"<<to_string(mode)<<endl;
	cout << "Number of variables - " <<num_Of_Variables <<endl;
	cout << "Number of Objective functions - " <<num_Of_Objective_Functions <<endl;
	cout << "Maximum deviation in converge limit " << errorlimit <<endl;
	cout << "Number of evaluations taken to converge - " <<evals <<endl;
	cout << "Total execution time : " << secs << "s" << endl;
	cout << "Variables values have been written to file Var_DTLZ"<<to_string(mode)<<"_"
			<<to_string(num_Of_Objective_Functions)<<".txt"<<endl;
	cout << "Objectives values have been written to file FUN_DTLZ"<<to_string(mode)<< "_"
			<<to_string(num_Of_Objective_Functions)<<".txt"<<endl;

	return 0;

}