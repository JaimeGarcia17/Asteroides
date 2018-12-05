#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <tgmath.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <chrono>
#include <omp.h>
#include <string.h>
#include <errno.h>


using namespace std;
using namespace std::chrono;
using clk = chrono::high_resolution_clock;


struct asteroides{
	double x,y,masa,vx,vy,fx,fy,ax,ay;
};

struct planetas{
	double x,y,masa;
};

struct rayo{
	double x,y,sup,inf;
};


/************************Comentarios**********************************
	El timer esta comentado, utilizo time en la linea de 
	comandos para comparar rendimiento, aun asi el timer es 
	mas preciso.

	El numero de hilos esta por defecto a 4, he dejado 
	comentado un pequeño trozo de codigo para que te pida
	el numero de hilos a ejecutar por pantalla, ya que,
	el enunciado pide que no se muestre nada por la misma, sólo
	habría que comentar el setnumthreads y descomentar las
	lineas siguientes, hay que tener en cuenta que el tiempo que se
	tarda en introducir el num de hilos tambien cuenta, por eso, 
	las pruebas las hemos realizado cambiando el valor en la sentencia
	set_num_threads.
*********************************************************************/



int main(int argc, char **argv) {


	/*Comienza timer*/

	/*auto t1 = clk::now();*/

	/*Constantes*/

	const double gravity = 6.674e-5;
	const double dt = 0.1;
	const double dmin = 2.0;
	const int width = 200;
	const int height = 200;
	const int ray_width = 4;
	const int m = 1000;
	const int sdm = 50;


	/*Comprobacion de parametros*/

	if (argc != 6 || atoi(argv[1]) < 0 || atoi(argv[2]) < 0 || atoi(argv[3]) < 0 || atof(argv[4]) < 0 || atof(argv[4]) > 200 || atoi(argv[5]) <= 0){
		cerr << "nasteroids-seq: Wrong arguments." <<endl <<"Correct use:" <<endl <<"nasteroids-seq num_asteroides num_iteraciones num_planetas pos_rayo semilla" <<endl;
		return -1;
	}

	int num_asteroides = atoi(argv[1]);
	int num_iteraciones = atoi(argv[2]);
	int num_planetas = atoi(argv[3]);
	double pos_rayo = atof(argv[4]);
	int semilla = atoi(argv[5]);



	/*Valores aleatorios*/

	default_random_engine re{semilla};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(width,std :: numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(height,std :: numeric_limits<double>::max())};
	normal_distribution<double> mdist{m, sdm};


	/*Declaracion de estructuras*/

	vector <asteroides> miasteroide;
	vector <planetas> miplaneta;
	rayo mirayo;


	/*Ficheros de salida*/
	ofstream fichero_inicial("init_conf.txt"); 
	ofstream fichero_final("out.txt"); 

	/*Parametros en fichero inicial*/

	fichero_inicial <<fixed <<setprecision(3) <<num_asteroides <<" " <<num_iteraciones <<" " <<num_planetas <<" " <<pos_rayo <<" " <<semilla <<endl;

	/*Colocacion asteroides*/

	for(int i = 0;i<num_asteroides; i++){
		miasteroide.push_back(asteroides());
		miasteroide[i].x = xdist(re);
		miasteroide[i].y = ydist(re);
		miasteroide[i].masa = mdist(re);
		fichero_inicial <<fixed <<setprecision(3) <<miasteroide[i].x <<" " <<miasteroide[i].y <<" " <<miasteroide[i].masa <<endl;
	}



	/*Colocacion planetas*/

	int contador=0;

	/*Las masas en teoría no se imprimen pero en su nasteroid-base si lo hace*/

	for(int i = 0; i<num_planetas;i++){

		miplaneta.push_back(planetas());

		if(contador==0){
			miplaneta[i].x = 0;
			miplaneta[i].y  = ydist(re);
			miplaneta[i].masa  = mdist(re)*10;
			contador++;
			contador=contador%4;
			fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<" " <<miplaneta[i].masa  <<endl;
			//fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<endl;
			continue;
		}
		if(contador==1){
			miplaneta[i].x  = xdist(re);
			miplaneta[i].y  = 0;
			miplaneta[i].masa  = mdist(re)*10;
			contador++;
			contador=contador%4;
			fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<" " <<miplaneta[i].masa  <<endl;
			//fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<endl;
			continue;	
		}
		if(contador==2){
			miplaneta[i].x  = width;
			miplaneta[i].y  = ydist(re);
			miplaneta[i].masa  = mdist(re)*10;
			contador++;
			contador=contador%4;
			fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<" " <<miplaneta[i].masa  <<endl;
			//fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<endl;
			continue;
		}
		if(contador==3){
			miplaneta[i].x  = xdist(re);
			miplaneta[i].y  = height;
			miplaneta[i].masa  = mdist(re)*10;
			contador++;
			contador=contador%4;
			fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<" " <<miplaneta[i].masa  <<endl;
			//fichero_inicial <<fixed <<setprecision(3) <<miplaneta[i].x  <<" " <<miplaneta[i].y  <<endl;
			continue;
		}
		
	}

	/*Colocacion rayo*/

	mirayo.x = 0.0;
	mirayo.y = pos_rayo;
	mirayo.sup = pos_rayo+(ray_width/2);
	mirayo.inf = pos_rayo-(ray_width/2);
	fichero_inicial <<fixed <<setprecision(3) <<mirayo.x <<" " <<mirayo.y <<endl;


	/*Variables*/

	double distanciaAst;
	double pendienteAst;
	double anguloAst;

	double distanciaPl;
	double pendientePl;
	double anguloPl;

	double fuerza;

	unsigned int i,j;

	/*Matrices dinámicas*/

	int fila_dinamica=num_asteroides;
	int columna_dinamica=num_asteroides;
	int columna_dinamica2=num_planetas;

	double **fuerzax=(double**)malloc(sizeof(double)*fila_dinamica);
	double **fuerzay=(double**)malloc(sizeof(double)*fila_dinamica);
	double **fuerzaxp=(double**)malloc(sizeof(double)*fila_dinamica);
	double **fuerzayp=(double**)malloc(sizeof(double)*fila_dinamica);

	for(int i=0;i<fila_dinamica;i++){
		fuerzax[i]=(double*)malloc(sizeof(double)*columna_dinamica);
		fuerzay[i]=(double*)malloc(sizeof(double)*columna_dinamica);
		fuerzaxp[i]=(double*)malloc(sizeof(double)*columna_dinamica2);
		fuerzayp[i]=(double*)malloc(sizeof(double)*columna_dinamica2);
	}


	/*Fijamos numero de hilos*/

	omp_set_num_threads(4);


	/*int hilos;
	cout <<"Por favor ingrese el número de hilos: ";
	cin >> hilos;
	omp_set_num_threads(hilos);*/

	/*Loop iteraciones*/

	for(int z=0;z<num_iteraciones;z++){

		/*Reseteo fuerzas de estructuras y matrices*/
		for(unsigned int i=0;i<miasteroide.size();i++){
			for(unsigned j=0;j<miasteroide.size();j++){		
				fuerzax[i][j]=0;
				fuerzay[i][j]=0;
			}
			miasteroide[i].fx=0;
			miasteroide[i].fy=0;
		}


		for(unsigned int i=0;i<miasteroide.size();i++){
			for(unsigned j=0;j<miplaneta.size();j++){		
				fuerzaxp[i][j]=0;
				fuerzayp[i][j]=0;
			}
		}

		
		/*Fuerzas gravitatorias*/

		/*Loop asteroide-asteroide*/

		for(i = 0; i < miasteroide.size(); i++){
			#pragma omp parallel for shared(i,miasteroide,fuerzax,fuerzay) private(j,fuerza,distanciaAst,pendienteAst,anguloAst)
			for(j = i+1;j < miasteroide.size(); j++){

				distanciaAst = sqrt(pow(miasteroide[i].x-miasteroide[j].x,2) + pow(miasteroide[i].y-miasteroide[j].y,2));

				if(distanciaAst>dmin){		

					pendienteAst = (miasteroide[i].y - miasteroide[j].y)/(miasteroide[i].x - miasteroide[j].x);

					if(pendienteAst >= 1 || pendienteAst <= -1){
						pendienteAst = pendienteAst - trunc(pendienteAst);
					} 

					anguloAst = atan(pendienteAst);

					fuerza=(gravity*miasteroide[i].masa*miasteroide[j].masa)/(pow(distanciaAst,2));

					if(fuerza>200){
						fuerza=200;
					}

					fuerzax[i][j]=fuerza*cos(anguloAst);
					fuerzay[i][j]=fuerza*sin(anguloAst);

					fuerzax[j][i]=-fuerzax[i][j];
					fuerzay[j][i]=-fuerzay[i][j];


				}

			}
		}



		/*Loop asteroide-planeta*/

		for(i = 0; i < miasteroide.size(); i++){
			#pragma omp parallel for shared(i,miasteroide,miplaneta,fuerzaxp,fuerzayp) private(j,fuerza,distanciaPl,pendientePl,anguloPl)
			for(j = 0;j < miplaneta.size(); j++){

				distanciaPl = sqrt(pow(miasteroide[i].x-miplaneta[j].x,2) + pow(miasteroide[i].y-miplaneta[j].y,2));
				
				pendientePl = (miasteroide[i].y - miplaneta[j].y)/(miasteroide[i].x - miplaneta[j].x);

				if(pendientePl >= 1 || pendientePl <= -1){
					pendientePl = pendientePl - trunc(pendientePl);
				} 

				anguloPl = atan(pendientePl);

				fuerza=(gravity*miasteroide[i].masa*miplaneta[j].masa)/(pow(distanciaPl,2));

				if(fuerza>200){
					fuerza=200;
				}

				fuerzaxp[i][j] = fuerza*cos(anguloPl);
				fuerzayp[i][j] = fuerza*sin(anguloPl);
				
			}

		}


		/*Acumulo fuerzas*/

		for(i=0;i<miasteroide.size();i++){
			for(j=0;j<miasteroide.size();j++){
				miasteroide[i].fx+=fuerzax[i][j];
				miasteroide[i].fy+=fuerzay[i][j];
			}
		}


		for(i=0;i<miasteroide.size();i++){
			for(j=0;j<miplaneta.size();j++){
				miasteroide[i].fx+=fuerzaxp[i][j];
				miasteroide[i].fy+=fuerzayp[i][j];
			}
		}


		/*Actualizo posiciones*/

		for(unsigned int i=0;i<miasteroide.size();i++){


			miasteroide[i].ax = miasteroide[i].fx/miasteroide[i].masa;
			miasteroide[i].ay = miasteroide[i].fy/miasteroide[i].masa;

			
			miasteroide[i].vx += (dt*miasteroide[i].ax);
			miasteroide[i].vy += (dt*miasteroide[i].ay);

			miasteroide[i].x += (dt*miasteroide[i].vx);
			miasteroide[i].y += (dt*miasteroide[i].vy);

			if(miasteroide[i].y <= 0){
				miasteroide[i].y = dmin;
				miasteroide[i].vy = -miasteroide[i].vy;
				continue;
			}

			if(miasteroide[i].y >= height){
				miasteroide[i].y = height - dmin; 
				miasteroide[i].vy = -miasteroide[i].vy;
				continue;
			}	

			if(miasteroide[i].x <= 0){
				miasteroide[i].x = dmin;
				miasteroide[i].vx = -miasteroide[i].vx; 
				continue;	
			}

			if(miasteroide[i].x >= width){
				miasteroide[i].x = height - dmin; 	
				miasteroide[i].vx = -miasteroide[i].vx; 
				continue;
			}



		}



		/*Rayo*/


		for(unsigned int i=0;i<miasteroide.size();i++){

			if(((miasteroide[i].y) <= mirayo.sup) && ((miasteroide[i].y) >= mirayo.inf)){
				miasteroide.erase(miasteroide.begin() +i);
			}

		}



	} /*Cierro loop iteraciones*/


	for(unsigned int i=0;i<miasteroide.size();i++){

		fichero_final <<fixed <<setprecision(3) <<miasteroide[i].x <<" " <<miasteroide[i].y <<" " <<miasteroide[i].vx <<" " <<miasteroide[i].vy <<" " <<miasteroide[i].masa <<endl;

	}


	/*Cierro ficheros*/

	fichero_inicial.close();
	fichero_final.close();

	/*Borro vectores*/

	miasteroide.clear();
	miplaneta.clear();

	/*Libero memoria de matrices*/

	for(int i=0;i<fila_dinamica;i++){
		free(fuerzax[i]);
		free(fuerzay[i]);
		free(fuerzaxp[i]);
		free(fuerzayp[i]);
	}

	free(fuerzax);
	free(fuerzay);
	free(fuerzaxp);
	free(fuerzayp);


	/*Apago timer*/

	/*
	auto t2 = clk::now();

	auto t3 = t2-t1;

	auto diff  = duration_cast<microseconds>(t3);

	cout <<diff.count() <<" ms" <<endl;*/

	return 0;
}


