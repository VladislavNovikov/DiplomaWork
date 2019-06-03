#include <iostream>
#include <cmath>
#include <fstream>
#include <windows.h>

void f(double a, double b, double c, double& RE, double& PE, double& RI) {
	RE = a;
	PE = b;
	RI = c;
}

using namespace std;

int main(){
	ofstream output;
	FILE* gnuplotPipe;
	gnuplotPipe = _popen("C:/gnuplot/bin/gnuplot.exe", "w");
	printf("Please press Enter");
	fflush(gnuplotPipe);

	const int n = 151;
	const int m = 151;
	double xmax, ymax, dx, dy, x[n], y[m];
	xmax = 10; dx = xmax/(n - 1);
	ymax = 1; dy = ymax/(m - 1);
	size_t i, j;

	for (i = 0; i < n; i++) { x[i] = dx * i;} for (j = 0; j < m; j++) { y[j] = dy * j;}

	// fprintf(gnuplotPipe, "set size ratio %lf\n", ymax/xmax);
	// fprintf(gnuplotPipe, "set xrange [%lf:%lf]\n", 0., xmax);
	// fprintf(gnuplotPipe, "set yrange [%lf:%lf]\n", 0., ymax);
	// fflush(gnuplotPipe);
	getchar();

	double **ksi_1 = new double*[n], **psi_1 = new double*[n], **theta1 = new double*[n], **u1 = new double*[n], **v1 = new double*[n];
	double **ksi_2 = new double*[n], **psi_2 = new double*[n], **theta2 = new double*[n], **u2 = new double*[n], **v2 = new double*[n];

	for (i = 0; i < n; i++) {
		ksi_1[i] = new double[m]; psi_1[i] = new double[m]; theta1[i] = new double[m]; u1[i] = new double[m]; v1[i] = new double[m];
		ksi_2[i] = new double[m]; psi_2[i] = new double[m]; theta2[i] = new double[m]; u2[i] = new double[m]; v2[i] = new double[m];
	}

	double ***ksi = new double**[2], ***psi = new double**[2], ***theta = new double**[2], ***u = new double**[2], ***v = new double**[2];  // ��� ��������
	ksi[0] = ksi_1; psi[0] = psi_1; theta[0] = theta1; u[0] = u1;	v[0] = v1;
	ksi[1] = ksi_2; psi[1] = psi_2; theta[1] = theta2; u[1] = u2;	v[1] = v2;

	double Th0 = 10;
	for (i = 0; i < n; i++) {
		//theta1[i][0] = 10.*exp(-pow((i*dx - dx * n/2), 2)/0.02);
    //theta2[i][0] = 10.*exp(-pow((i*dx - dx * n/2), 2)/0.02);


		theta1[i][0] = Th0 + 0*sin(10*3.141592*i/n) + 0.1*(double)rand()/RAND_MAX;
		//theta1[i][0] = Th0 + sin(8*3.141592*i/n);
		theta2[i][0] = theta1[i][0];

	}

	double u0 = 1, L = 200, dt = 0.001, t0 = L/u0; // ��������� �������

	double alfa, beta, nu, Re, Pe, Ri, Ra, xi, g;

	beta = 0.0001; xi = 0.02; nu = 0.02; g = 10;
	Re = u0 * L / nu; Pe = u0 * L / xi; Ra = g * beta * Th0 * L*L*L / (nu*xi); Ri = Ra / (Re * Pe);

//	f(100, 50, 10, Re, Pe, Ri);
	//f(1000, 1000, 10, Re, Pe, Ri);
	//f(10000, 10000, 10, Re, Pe, Ri);


	int s = 0, s_max = 1000000, K, s_image = 100;
	double eps = 0.01, max;

	while (s < s_max) {
    // явная схема ksi
		for (i = 1; i < n - 1; i++) {
			for (j = 1; j < m - 1; j++) {
					if (u[s%2][i][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][i][j] >= 0) { beta = 0; } else { beta = 1; }

					ksi[(s+1)%2][i][j] = ksi[s%2][i][j]
						- dt/dx * u[s%2][i][j] * ((1 - alfa)*(ksi[s%2][i][j] - ksi[s%2][i-1][j]) + alfa*(ksi[s%2][i+1][j] - ksi[s%2][i][j]))
						- dt/dy * v[s%2][i][j] * ((1 - beta)*(ksi[s%2][i][j] - ksi[s%2][i][j-1]) + beta*(ksi[s%2][i][j+1] - ksi[s%2][i][j]))
							+ dt/Re*((ksi[s%2][i+1][j] - 2 * ksi[s%2][i][j] + ksi[s%2][i-1][j])/dx/dx
							+ (ksi[s%2][i][j+1] - 2 * ksi[s%2][i][j] + ksi[s%2][i][j-1])/dy/dy)
								- dt * Ri/(2*dx) * (theta[s%2][i+1][j] - theta[s%2][i-1][j]); }}

		/////// ГРАНИЦА ksi слева и справа (периодические граничные условия).
		for (j = 1; j < m - 1; j++) {
				if (u[s%2][0][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][0][j] >= 0) { beta = 0; } else { beta = 1; }

				ksi[(s+1)%2][0][j] = ksi[s%2][0][j]
					- dt/dx * u[s%2][0][j] * ((1 - alfa)*(ksi[s%2][0][j] - ksi[s%2][n-2][j]) + alfa*(ksi[s%2][1][j] - ksi[s%2][0][j]))
					- dt/dy * v[s%2][0][j] * ((1 - beta)*(ksi[s%2][0][j] - ksi[s%2][0][j-1]) + beta*(ksi[s%2][0][j+1] - ksi[s%2][0][j]))
						+ dt/Re*((ksi[s%2][1][j] - 2 * ksi[s%2][0][j] + ksi[s%2][n-2][j])/dx/dx
						+ (ksi[s%2][0][j+1] - 2 * ksi[s%2][0][j] + ksi[s%2][0][j-1])/dy/dy)
							- dt * Ri/(2*dx) * (theta[s%2][1][j] - theta[s%2][n-2][j]);
				ksi[(s+1)%2][n-1][j] = ksi[(s+1)%2][0][j]; // справа = слева
		}

    int K = 0;
		do { // Пуассон
			max = 0;
			// для произвольной формы. внутри
			for (i = 1; i < n - 1; i++){
				for (j = 1; j < m - 1; j++){
						psi[(s+1) % 2][i][j] = (dy*dy*(psi[s%2][i+1][j] + psi[(s+1)%2][i-1][j]) + dx * dx*(psi[s%2][i][j+1] + psi[(s+1)%2][i][j-1])
							- ksi[s%2][i][j] * dx*dx*dy*dy)
							/ (2. * (dy * dy + dx * dx));
			}}
			////// ГРАНИЦА psi слева и справа (периодические граничные условия)
			for (j = 1; j < m - 1; j++){
						psi[(s+1) % 2][0][j] = (dy*dy*(psi[s%2][1][j] + psi[(s+1) % 2][n-2][j]) + dx * dx*(psi[s%2][0][j+1] + psi[(s+1)% 2][0][j-1])
							- ksi[s%2][0][j] * dx*dx*dy*dy)
							/ (2. * (dy*dy + dx * dx));

						psi[(s+1) % 2][n-1][j] = psi[(s+1) % 2][0][j]; // справа = слева
			}
			for ( i = 1; i < n - 1; i++){
				for ( j = 1; j < m - 1; j++){
					max += abs(psi[(s+1)%2][i][j] - psi[s%2][i][j]);
					psi[s%2][i][j] = psi[(s+1)%2][i][j]; // перезапись, а то не пуассон не сойдётся
			}}

			K++;
		}
		while (max > eps);

    // Скорости
       for ( i = 1; i < n - 1; i++){
         for ( j = 1; j < m - 1; j++){
           // производная по 2м точкам
           u[(s+1)%2][i][j] = (psi[(s+1)%2][i][j + 1] - psi[(s+1)%2][i][j - 1]) / 2. / dy;
           v[(s+1)%2][i][j] = -(psi[(s+1)%2][i + 1][j] - psi[(s+1)%2][i - 1][j]) / 2. / dx;
         }}

         // скорости (периодические граничные условия)
       for (j = 1; j < m - 1; j++){
         // производная по 2м соседним точкам
         u[(s+1)%2][0][j] = (psi[(s+1)%2][0][j+1] - psi[(s+1)%2][0][j-1])/2./dy;
         u[(s+1)%2][n-1][j] = u[(s+1)%2][0][j]; // справа = слева

         v[(s+1)%2][0][j] = -(psi[(s+1)%2][1][j] - psi[(s+1)%2][n-2][j])/2./dx;
         v[(s+1)%2][n-1][j] = v[(s+1)%2][0][j]; // справа = слева

       // или
       // 	u[(s+1)%2][n - 1][j] = u[(s+1)%2][n - 2][j];
       // 	u[(s+1)%2][0][j] = u[(s+1)%2][n-1][j];
       // 	v[(s+1)%2][n - 1][j] = v[(s+1)%2][n - 2][j];
       // 	v[(s+1)%2][0][j] = v[(s+1)%2][n-1][j];
      }

      ///////////////////// РАСЧЁТ температуры theta
   	 for (i = 1; i < n - 1; i++) {
   		 for (j = 1; j < m - 1; j++){
   				 if (u[s%2][i][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][i][j] >= 0) { beta = 0; } else { beta = 1; }

   				 theta[(s+1)%2][i][j] = theta[s%2][i][j]
   					 - dt/dx * u[s%2][i][j] * ((1 - alfa)*(theta[s%2][i][j] - theta[s%2][i-1][j]) + alfa*(theta[s%2][i+1][j] - theta[s%2][i][j]))
   					 - dt/dy * v[s%2][i][j] * ((1 - beta)*(theta[s%2][i][j] - theta[s%2][i][j-1]) + beta*(theta[s%2][i][j+1] - theta[s%2][i][j]))
   						 +(1/Pe)*dt * ((theta[s%2][i+1][j] - 2 * theta[s%2][i][j] + theta[s%2][i-1][j]) /dx/dx
   							 + (theta[s%2][i][j+1] - 2 * theta[s%2][i][j] + theta[s%2][i][j-1]) /dy/dy);
   		 }
   	 }
   	 //// границы по theta слева и справа (периодические граничные условия)
   	 for (j = 1; j < m-1; j++) {
   		 if (u[s%2][0][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][0][j] >= 0) { beta = 0; } else { beta = 1; }

   		 theta[(s+1)%2][0][j] = theta[s%2][0][j]
   			 - dt/dx * u[s%2][0][j] * ((1 - alfa)*(theta[s%2][0][j] - theta[s%2][n-2][j]) + alfa*(theta[s%2][1][j] - theta[s%2][0][j]))
   			 - dt/dy * v[s%2][0][j] * ((1 - beta)*(theta[s%2][0][j] - theta[s%2][0][j-1]) + beta*(theta[s%2][0][j+1] - theta[s%2][0][j]))
   				 +(1/Pe)*dt * ((theta[s%2][1][j] - 2 * theta[s%2][0][j] + theta[s%2][n-2][j]) /dx/dx
   					 + (theta[s%2][0][j+1] - 2 * theta[s%2][0][j] + theta[s%2][0][j-1]) /dy/dy);
   		 theta[(s+1)%2][n-1][j] = theta[(s+1)%2][0][j]; // справа = слева
   	 }

		//////////////////////////////////////// ������� (psi, ksi) ������� ������� � ������. psi = 0;
		for (i = 0; i < n; i++) {
			 ksi[(s+1)%2][i][m - 1] = 2 * (psi[(s+1)%2][i][m - 2] - psi[(s+1)%2][i][m - 1])/dy/dy;
			 ksi[(s+1)%2][i][0] = 2 * (psi[(s+1)%2][i][1] - psi[(s+1)%2][i][0])/dy/dy;
			//
			// u[(s+1)%2][i][m-1] = (psi[(s+1)%2][i][m-1] - psi[(s+1)%2][i][m-2])/2./dy;
			// u[(s+1)%2][i][0] = (psi[(s+1)%2][i][1] - psi[(s+1)%2][i][0])/2./dy;

		/////////////////////////////////////// ������� �� theta ������ � �����(������)
			theta[(s+1)%2][i][m - 1] = theta[(s+1)%2][i][m - 2] * (1. - 0.2);
			//theta[(s+1)%2][i][0] = Th0 + 0*sin(10*3.141592*i/n) + 0.1*(double)rand()/RAND_MAX;
			theta[(s+1)%2][i][0] = Th0 + 0.1*(double)rand()/RAND_MAX;
		}

		// запись в файл
		if (s % s_image == 0) {
			fprintf(gnuplotPipe, "plot '-' using 1:2:3 with image \n");
			//output.open("FILE.txt");
			for (i = 0; i < n; i++){
				for (j = 0; j < m; j++){
				//	output << x[i] << "  " << y[j] << "  " << theta[(s+1)%2][i][j] << endl;

					fprintf(gnuplotPipe, "%lf, %lf, %lf\n", x[i],y[j], theta[(s+1)%2][i][j]	);
					//fprintf(gnuplotPipe, "%lf, %lf, %lf\n", x[i],y[j], v[(s+1)%2][i][j]	);
				//	fprintf(gnuplotPipe, "%lf, %lf, %lf\n", x[i],y[j], sqrt(pow(u[(s+1)%2][i][j],2) + pow(v[(s+1)%2][i][j], 2))	);

					//std::cout << theta[(s+1)%2][i][j] << '\n';
					//output << x[i] << "  " << y[j] << "  " << v[(s+1)%2][i][j] << endl;
					//output << x[i] << "  " << y[j] << "  " << sqrt(pow(u[(s+1)%2][i][j],2) + pow(v[(s+1)%2][i][j], 2)) << endl;
					//output << x[i] << "  " << y[j] << "  " << v[(s+1)%2][i][j]<<endl;
					//output << x[i] << "  " << y[j] << "  " << u[(s+1)%2][i][j] << " " << v[(s+1)%2][i][j] << " "
						//<< sqrt(u[(s+1)%2][i][j] * u[(s+1)%2][i][j] + v[(s+1)%2][i][j] * v[(s+1)%2][i][j]) << endl;
				}
			}

			//for (i = 0; i < n; i++){	output << i * dx << "  " << v[(s+1)%2][i][Gran] << endl; } // одномерный
			//fprintf(gnuplotPipe, "plot 'FILE.txt' with lines \n");

			//output.close();
			system("cls");
			cout << "Lx=" << xmax*L << " Ly=" << ymax*L <<endl;
			cout << "Ra=" << Ra <<endl;
			cout << "Re=" << Re << " Pe=" << Pe << " Ri=" << Ri << endl;
			cout << "K=" << K <<" step=" << s << " time=" << dt * (s+1) * t0 / 60<< " min" << endl;

			//fprintf(gnuplotPipe, "plot 'FILE.txt' using 1:2:3 with image \n");
			//fprintf(gnuplotPipe, "plot 'FILE.txt' with vectors lc palette \n");
			//fprintf(gnuplotPipe, "plot '%s' u 1:2:($3*0.1):($4*0.1):5 with vectors head size 0.,0,0 filled lc palette\n", "FILE.txt");
			fprintf(gnuplotPipe, "e\n");
			fflush(gnuplotPipe);

			//getchar();
		}
		s++;
	}
	cout << endl;
	printf("press enter to close gnuplot\n");
	getchar();
	fprintf(gnuplotPipe, "exit\n");
	return 0;
	}
