#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "seal.h"
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;
using namespace seal;
/*float calcul_Erreur(float theta_0, float theta_1, float theta_2, float theta_3, float theta_4, float theta_5, float theta_6, float theta_7, float theta_8,vector<float> & V1,vector<float> & V2,vector<float> & V3,vector<float> & V4,vector<float> & V5,vector<float> & V6,vector<float> & V7,vector<float> & V8,vector<float> & classe, int N){
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
			}*/
std::vector<float> stochastic_gradient_descent_seal(vector<float> & V1,vector<float> & V2,vector<float> & V3,vector<float> & V4,vector<float> & V5,vector<float> & V6,vector<float> & V7,vector<float> & V8,vector<float> & classe, int N){
  //vector <float> theta{1,1,1,1,1,1,1,1,1};
  float theta_0 = 1,theta_1 = 1, theta_2 = 1, theta_3 = 1, theta_4 = 1, theta_5 = 1, theta_6 = 1, theta_7 = 1, theta_8 = 1;
         float Erreur;
         float learningRate = 0.001;
         int NumberOfIter=0;
         int i = 0;
         int iter=0;
	 EncryptionParameters parms;
	 parms.set_poly_modulus("1x^2048 + 1");
	 parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(2048));
	 parms.set_plain_modulus(1 << 4);
	 parms.set_decomposition_bit_count(16);
	 parms.validate();
	 KeyGenerator generator(parms);
	 generator.generate(1);
	 Ciphertext public_key = generator.public_key();
	 Plaintext secret_key = generator.secret_key();
	 EvaluationKeys evaluation_keys = generator.evaluation_keys();
	 FractionalEncoder encoder(parms.plain_modulus(), parms.poly_modulus(), 64, 32, 3);
	 Encryptor encryptor(parms, public_key);
	 Decryptor decryptor(parms, secret_key);
	 Evaluator evaluator(parms, evaluation_keys);
	 Plaintext encodedLearningRate = encoder.encode(learningRate);
	 Plaintext encodedTheta_0 = encoder.encode(theta_0);
	 Plaintext encodedTheta_1 = encoder.encode(theta_1);
	 Plaintext encodedTheta_2 = encoder.encode(theta_2);
	 Plaintext encodedTheta_3 = encoder.encode(theta_3);
	 Plaintext encodedTheta_4 = encoder.encode(theta_4);
	 Plaintext encodedTheta_5= encoder.encode(theta_5);
	 Plaintext encodedTheta_6 = encoder.encode(theta_6);
	 Plaintext encodedTheta_7 = encoder.encode(theta_7);
	 Plaintext encodedTheta_8 = encoder.encode(theta_8);
	 
	 vector<Ciphertext> encrypted_V_1;
	 vector<Ciphertext> encrypted_V_2;
	 vector<Ciphertext> encrypted_V_3;
	 vector<Ciphertext> encrypted_V_4;
	 vector<Ciphertext> encrypted_V_5;
	 vector<Ciphertext> encrypted_V_6;
	 vector<Ciphertext> encrypted_V_7;
	 vector<Ciphertext> encrypted_V_8;
	 vector<Ciphertext> encrypted_classe;
	 vector<Plaintext> encodedTheta;
	 vector<Ciphertext> somme_mult;
	 vector<Ciphertext> tay;
	 //Plaintext encodedTheta_0 = encoder.encode(theta_0);
	 /*for (int i = 0; i < 9; ++i){
	   Plaintext encoded_theta = encoder.encode(theta[i]);
	   encodedTheta.emplace_back(encoded_theta);
	   }*/
	 
	 for (int i = 0; i < N; ++i){
	   Plaintext encoded_V_1 = encoder.encode(V1[i]);
	   Plaintext encoded_V_2 = encoder.encode(V2[i]);
	   Plaintext encoded_V_3 = encoder.encode(V3[i]);
	   Plaintext encoded_V_4 = encoder.encode(V4[i]);
	   Plaintext encoded_V_5 = encoder.encode(V5[i]);
	   Plaintext encoded_V_6 = encoder.encode(V6[i]);
	   Plaintext encoded_V_7 = encoder.encode(V7[i]);
	   Plaintext encoded_V_8 = encoder.encode(V8[i]);
	   Plaintext encoded_classe = encoder.encode(classe[i]);
	   encrypted_V_1.emplace_back(encryptor.encrypt(encoded_V_1));
	   encrypted_V_2.emplace_back(encryptor.encrypt(encoded_V_2));
	   encrypted_V_3.emplace_back(encryptor.encrypt(encoded_V_3));
	   encrypted_V_4.emplace_back(encryptor.encrypt(encoded_V_4));
	   encrypted_V_5.emplace_back(encryptor.encrypt(encoded_V_5));
	   encrypted_V_6.emplace_back(encryptor.encrypt(encoded_V_6));
	   encrypted_V_7.emplace_back(encryptor.encrypt(encoded_V_7));
	   encrypted_V_8.emplace_back(encryptor.encrypt(encoded_V_8));
	   encrypted_classe.emplace_back(encryptor.encrypt(encoded_classe));
	 }
	
	 while (iter < 1){
         	 for (i=0 ; i < 2 ; i++){
		   cout<<"("<<i+1<<") tour de boucle"<<endl;
		   encodedTheta_0 = encoder.encode(theta_0);
		   cout<<"theta_0 = "<<theta_0<<endl;
		   
		   somme_mult.clear();
		   /*begin sigmoid function*/
		   Ciphertext mult_theta1 = evaluator.multiply_plain(encrypted_V_1[i],encodedTheta_1);
		   somme_mult.emplace_back(mult_theta1);
		   Ciphertext mult_theta2 = evaluator.multiply_plain(encrypted_V_2[i],encodedTheta_2);
		   somme_mult.emplace_back(mult_theta2);
		   Ciphertext mult_theta3 = evaluator.multiply_plain(encrypted_V_3[i],encodedTheta_3);
		   somme_mult.emplace_back(mult_theta3);
		   Ciphertext mult_theta4 = evaluator.multiply_plain(encrypted_V_4[i],encodedTheta_4);
		   somme_mult.emplace_back(mult_theta4);
		   Ciphertext mult_theta5 = evaluator.multiply_plain(encrypted_V_5[i],encodedTheta_5);
		   somme_mult.emplace_back(mult_theta5);
		   Ciphertext mult_theta6 = evaluator.multiply_plain(encrypted_V_6[i],encodedTheta_6);
		   somme_mult.emplace_back(mult_theta6);
		   Ciphertext mult_theta7 = evaluator.multiply_plain(encrypted_V_7[i],encodedTheta_7);
		   somme_mult.emplace_back(mult_theta7);
		   Ciphertext mult_theta8 = evaluator.multiply_plain(encrypted_V_8[i],encodedTheta_8);
		   somme_mult.emplace_back(mult_theta8);
		   Ciphertext somme = evaluator.add_many(somme_mult);
		   Plaintext test1 = decryptor.decrypt(somme);
		   float     test_1 = encoder.decode(test1);
		   cout<<"-----------------------------------------------------------"<<endl;
		   cout<<"test_1 = "<<test_1<<endl;
		   somme_mult.clear();
		   Ciphertext somme_theta0 = evaluator.add_plain(somme, encodedTheta_0);
		   Plaintext test2 = decryptor.decrypt(somme_theta0);
		   float     test_2 = encoder.decode(test2);
		   cout<<"test_2 = "<<test_2<<endl;
		   Ciphertext somme_negative = evaluator.negate(somme_theta0);
		   Plaintext test3 = decryptor.decrypt(somme_negative);
		   float test_3 = encoder.decode(test3);
		   cout <<"test_3 = "<<test_3<<endl;
		   vector<Plaintext> taylor_coeffs {encoder.encode(1.0/4),encoder.encode(-1.0/48),encoder.encode(1.0/480),encoder.encode(-17.0/80640)};
		   
		   Plaintext taylor_constant = encoder.encode(1.0/2);
		   Ciphertext result1 = evaluator.multiply_plain(somme_negative, taylor_coeffs[0]);
		   Plaintext test4 = decryptor.decrypt(result1);
		   float test_4 = encoder.decode(test4);
		   cout <<"test_4 = "<<test_4<<endl;
		   tay.emplace_back(result1);
		   Ciphertext result2 = evaluator.multiply_plain(evaluator.exponentiate(somme_negative,1), taylor_coeffs[1]);
		   Plaintext test5 = decryptor.decrypt(result2);
		   float test_5= encoder.decode(test5);
		   cout <<"test_5 = "<<test_5<<endl;
		   
		   tay.emplace_back(result2);
		   Ciphertext result3 = evaluator.multiply_plain(evaluator.exponentiate(somme_negative, 5), taylor_coeffs[2]);
		   Plaintext test6 = decryptor.decrypt(result3);
		   float test_6= encoder.decode(test6);
		   cout <<"test_6 = "<<test_6<<endl;
		   
		   tay.emplace_back(result3);
		   
		   Ciphertext result4 = evaluator.multiply_plain(evaluator.exponentiate(somme_negative,7), taylor_coeffs[3]);
		   Plaintext test7 = decryptor.decrypt(result4);
		   float test_7= encoder.decode(test7);
		   cout <<"test_7 = "<<test_7<<endl;
		   
		   tay.emplace_back(result4);
		   Ciphertext somme_tay = evaluator.add_many(tay);
		   Ciphertext resul = evaluator.add_plain(somme_tay, taylor_constant);
		   /* end sigmoid function*/
		   /*
		   Ciphertext diff = evaluator.sub(result, encrypted_classe[i]);
		   Plaintext inter1 = decryptor.decrypt(diff);
		   float inter_1 = encoder.decode(inter1);
		  
		   Ciphertext mult = evaluator.multiply_plain(diff, encodedLearningRate);
		   Ciphertext inter4 = evaluator.sub_plain(mult,  encodedTheta_0);
		   Ciphertext encryptedtheta_0 = evaluator.negate(inter4);
		   Plaintext encodedthe = decryptor.decrypt(encryptedtheta_0);
		   theta_0 = encoder.decode(encodedthe);
		   cout<<"theta_0 = "<<theta_0<<endl;*/
		   
		   /*  theta[0]=theta[0] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i]);*/
		     /* theta[1]=theta[1] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V1[i];
                     theta[2]=theta[2] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V2[i];
                     theta[3]=theta[3] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V3[i];
                     theta[4]=theta[4] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V4[i];
                     theta[5]=theta[5] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V5[i];
                     theta[6]=theta[6] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V6[i];
                     theta[7]=theta[7] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V7[i];
                     theta[8]=theta[8] - learningRate* ((1/(1+exp(-(theta[0] + theta[1]*V1[i] + theta[2]*V2[i]+theta[3]*V3[i]+theta[4]*V4[i]+theta[5]*V5[i]+theta                     [6]*V6[i]+theta[7]*V7[i]+theta[8]*V8[i]))))-classe[i])*V8[i];*/
            }
								  /* Erreur = Erreur = calcul_Erreur(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6],theta[7],theta[8], V1, V2, V3, V4, V5, V6, V7, V8, classe, N);
								     monFlux << iter << " , " << Erreur << " ; " << endl;*/
        iter = iter + 1;
        }
							    // return  theta;
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
        theta = stochastic_gradient_descent_seal(V1_F, V2_F, V3_F, V4_F, V5_F, V6_F, V7_F,V8_F, classe_F,V1_F.size());
	/*for(int i=V1_F.size()-201; i<V1_F.size(); i++){
	         temp = predict(theta,V1_F[i], V2_F[i], V3_F[i], V4_F[i], V5_F[i], V6_F[i], V7_F[i],V8_F[i]);
		 if(temp == classe_F[i])    t = t+1;
		 
	}
	cout<< "accuracy = "  << (t/201)*100 << " % " << endl;*/
	return 1;
}    
    
