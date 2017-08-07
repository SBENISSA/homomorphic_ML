#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;
float calcul_Erreur(float theta_0, float theta_1, float theta_2, float theta_3, float theta_4, float theta_5, float theta_6, float theta_7, float theta_8,vector<float> & V1,vector<float> & V2,vector<float> & V3,vector<float> & V4,vector<float> & V5,vector<float> & V6,vector<float> & V7,vector<float> & V8,vector<float> & classe, int N){
        float Erreur = 0;
        float a, b, c, d, e, f, g, h, k;
        for (int i=0; i < N ; i++){
              a = V1[i];
              b = V2[i];
              c = V3[i];
              d = V4[i];
              e = V5[i];
              f = V6[i];
              g = V7[i];
              h = V8[i];
              k = classe[i];
              Erreur += k*log (1/(1+exp(-(theta_0 + theta_1*a + theta_2 * b + theta_3*c,theta_4*d+ theta_5*e + theta_6 * f + theta_7*g,theta_8*h)))) + (1-k)*log               (1 - 1/(1+exp(-(theta_0 + theta_1*a + theta_2 * b + theta_3*c,theta_4*d+ theta_5*e + theta_6 * f + theta_7*g,theta_8*h)))); 
            }
         return -((Erreur/float(N)));
}
bool predict(vector<float> & theta, float X_1, float X_2, float X_3, float X_4, float X_5, float X_6, float X_7, float X_8){
          if ((1/(1+exp(-(theta[0]+theta[1]*X_1+theta[2]*X_2+theta[3]*X_3+theta[4]*X_4+theta[5]*X_5+theta[6]*X_6+theta[7]*X_7+theta[8]*X_8)))) > 0.5)
                     {
                        return 1;
                      }
          else 
                        return 0;
}
std::vector<float> stochastic_gradient_descent(vector<float> & V1,vector<float> & V2,vector<float> & V3,vector<float> & V4,vector<float> & V5,vector<float> & V6,vector<float> & V7,vector<float> & V8,vector<float> & classe, int N){
  vector <float> theta{1,1,1,1,1,1,1,1,1};
         float Erreur;
         float learningRate = 0.001;
	 // int NumberOfIter=0;
         int i = 0;
         int iter=0;
         string const nomFichier("erreur.csv");
         ofstream monFlux(nomFichier.c_str());
         while (iter < 1){
	   cout<<"------------------------------------------------------------------------------"<<endl;
         for (i=0 ; i<2 ; i++){
	   cout<<"------------------------------------------------------"<<endl;
	     float inter1 =theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta[6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i];
	     cout<<"inter1 = "<<inter1<<endl;
	     float inter2 =theta[0] + theta[1]*V1[i]+theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta[6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i];
	     cout<<"inter2 = "<<inter2<<endl;
	     float inter3 = -(theta[0]+ theta[1]*V1[i]+theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta[6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]);
	     cout <<"inter3 = "<<inter3<<endl;
	     float inter4 =  (float(1)/4)*inter3 ;
	     cout <<"inter4 = "<<inter4<<endl;
	     float inter5 =  (float(1)/48)*pow(inter3,3) ;
	     cout <<"inter5 = "<<inter5<<endl;
	     float inter6 =  (float(1/480))*pow(inter3,5) ;
	     cout <<"inter6 = "<<inter6<<endl;
	     float inter7 =  (float(17)/80640)*pow(inter3,3) ;
	     cout <<"inter7 = "<<inter7<<endl;
	    float inter5 = (1/(1+exp(-(theta[0]+theta[1]*V1[i]+theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta \
				       [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))));
	   cout <<"inter5 = "<<inter5<<endl;
	   float inter6 = ((1/(1+exp(-(theta[0]+theta[1]*V1[i]+theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta          \
				       [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i]);
	   cout<<"inter6 = "<<inter6<<endl;
	   float inter7 = learningRate* ((1/(1+exp(-(theta[0]+theta[1]*V1[i]+theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta          \
						     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i]);
	   cout<<"inter7 = "<<inter7<<endl;
	   theta[0]=theta[0] - learningRate* ((1/(1+exp(-(theta[0]+theta[1]*V1[i]+theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i]);
	   cout <<"theta[0] = "<<theta[0]<<endl;
	 }
		     /*theta[1]=theta[1] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V1[i];
                     theta[2]=theta[2] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V2[i];
                     theta[3]=theta[3] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V3[i];
                     theta[4]=theta[4] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V4[i];
                     theta[5]=theta[5] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V5[i];
                     theta[6]=theta[6] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V6[i];
                     theta[7]=theta[7] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V7[i];
                     theta[8]=theta[8] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V8[i];
            }
        Erreur = Erreur = calcul_Erreur(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6],theta[7],theta[8], V1, V2, V3, V4, V5, V6, V7, V8, classe, N);
        monFlux << iter << " , " << Erreur << " ; " << endl;*/
	 iter = iter + 1;
        }
        return  theta;
}
int main(){
	 vector <string> ensemble ;
         vector <string> V1_S;
         vector <string> V2_S;
         vector <string> V3_S;
         vector <string> V4_S;
	 vector <string> V5_S;
	 vector <string> V6_S;
	 vector <string> V7_S;
	 vector <string> V8_S;
	 vector <string> classe_S;
	 vector <float> V1_F;
	 vector <float> V2_F;
	 vector <float> V3_F;
	 vector <float> V4_F;
	 vector <float> V5_F;
	 vector <float> V6_F;
	 vector <float> V7_F;
	 vector <float> V8_F;
	 vector <float> classe_F;
	 vector <float> theta(9);
	 bool temp;
	 float t = 0;     
	 ifstream csv("pima.csv");
	 string line;
	 if (csv.is_open()) {
	        for (int i = 0; csv.good(); i++){
		            getline(csv, line, ',');
                            ensemble.push_back(line);
		}
	 }
	 else {
	      cout << "Unable to open file";
	 }
	 csv.close();
	 for (int i = 0; i != ensemble.size(); i++){
                       if ((i%9) == 0 )       V1_S.push_back(ensemble[i]);
		       else if ((i%9) == 1 )  V2_S.push_back(ensemble[i]);
		       else if ((i%9) == 2 )  V3_S.push_back(ensemble[i]);   
		       else if((i%9) == 3)    V4_S.push_back(ensemble[i]);
		       else if((i%9) == 4)    V5_S.push_back(ensemble[i]);
		       else if((i%9) == 5)    V6_S.push_back(ensemble[i]);
		       else if((i%9) == 6)    V7_S.push_back(ensemble[i]);
		       else if((i%9) == 7)    V8_S.push_back(ensemble[i]);
		       else   classe_S.push_back(ensemble[i]);
	 }
	 for (int i = 0; i != V1_S.size(); i++){
	              float a, b, c, d, e, f, g, h, k;
		      a = strtof(V1_S[i].c_str(),0);
		      b = strtof( V2_S[i].c_str(),0);
		      c = strtof( V3_S[i].c_str(),0);
		      d = strtof( V4_S[i].c_str(),0);
		      e = strtof( V5_S[i].c_str(),0);
		      f = strtof( V6_S[i].c_str(),0);
		      g = strtof( V7_S[i].c_str(),0);
		      h = strtof( V8_S[i].c_str(),0);
		      k = strtof( classe_S[i].c_str(),0);
		      V1_F.push_back(a);
		      V2_F.push_back(b);
		      V3_F.push_back(c);
		      V4_F.push_back(d);
		      V5_F.push_back(e);
		      V6_F.push_back(f);
		      V7_F.push_back(g);
		      V8_F.push_back(h);
		      classe_F.push_back(k);
	 }
	printf("stochastic gradient descent\n");
	theta = stochastic_gradient_descent(V1_F, V2_F, V3_F, V4_F, V5_F, V6_F, V7_F,V8_F, classe_F,V1_F.size());
	/*for(int i=V1_F.size()-201; i<V1_F.size(); i++){
	         temp = predict(theta,V1_F[i], V2_F[i], V3_F[i], V4_F[i], V5_F[i], V6_F[i], V7_F[i],V8_F[i]);
		 if(temp == classe_F[i])    t = t+1;
		 
	}
	cout<< "accuracy = "  << (t/201)*100 << " % " << endl;*/
	return 1;
}    
    
