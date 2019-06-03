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

	const int n = 241;
 	const int m = 81;
	int Gran = m/2.;
	double xmax, ymax, dx, dy, x[n], y[m];
	xmax = 3; dx = xmax/(n - 1);
	ymax = 1; dy = ymax/(m - 1);
	size_t i, j;

	for (i = 0; i < n; i++) { x[i] = dx * i;} for (j = 0; j < m; j++) { y[j] = dy * j;}

	fprintf(gnuplotPipe, "set size ratio %lf\n", ymax/xmax);
	fprintf(gnuplotPipe, "set xrange [%lf:%lf]\n", 0., xmax);
	fprintf(gnuplotPipe, "set yrange [%lf:%lf]\n", 0., ymax);
	fflush(gnuplotPipe);
	getchar();

	double **ksi_1 = new double*[n], **psi_1 = new double*[n], **theta1 = new double*[n], **u1 = new double*[n], **v1 = new double*[n];
  double **ksi_2 = new double*[n], **psi_2 = new double*[n], **theta2 = new double*[n], **u2 = new double*[n], **v2 = new double*[n];

	for (i = 0; i < n; i++) {
		ksi_1[i] = new double[m]; psi_1[i] = new double[m]; theta1[i] = new double[m]; u1[i] = new double[m]; v1[i] = new double[m];
		ksi_2[i] = new double[m]; psi_2[i] = new double[m]; theta2[i] = new double[m]; u2[i] = new double[m]; v2[i] = new double[m];
	}

	// ��� ��������
	double ***ksi = new double**[2], ***psi = new double**[2], ***theta = new double**[2], ***u = new double**[2], ***v = new double**[2];
	ksi[0] = ksi_1; psi[0] = psi_1; theta[0] = theta1; u[0] = u1;	v[0] = v1;
	ksi[1] = ksi_2; psi[1] = psi_2; theta[1] = theta2; u[1] = u2;	v[1] = v2;

	double T_1, T_2, Th0;
	Th0 = 10;	T_1 = 0;	T_2 = Th0;

	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++) {
			if (j <= Gran) {
					theta1[i][j] = T_2;	theta2[i][j] = T_2;
					v1[i][j] = 0.1*abs(powf(j/(double)Gran,4))*sin(10*3.141592*i/n); v2[i][j] = v1[i][j];
			 }
			else {
				theta1[i][j] = T_1;	theta2[i][j] = T_1;
			  v1[i][j] = v1[i][2*Gran - j];
			}
	}}
	for (i = 0; i < n ; i++){
		//theta1[i][Gran] = T_2 + 0.1*sin(10*3.141592*i/n);
	}

	double u0 = 1;
	for (i = 0; i < n; i++) {
		u1[i][m-1] = u0; u2[i][m-1] = u1[i][m-1];
		u1[i][0] = -u0; u2[i][0] = u1[i][0];
		}

  // ������
	for (i = 1; i < n - 1; i++) {
		for (j = 1; j < m - 1; j++) {
			ksi_1[i][j] = (u1[i][j+1] - u1[i][j-1]) / (2*dy) - (v1[i+1][j] - v1[i-1][j]) / (2*dx);
			ksi_2[i][j] = ksi_1[i][j];
		}}
	// ������ � �����
	for (j = 1; j < m - 1; j++) {
			ksi_1[0][j] = (u1[0][j+1] - u1[0][j-1]) / (2*dy) - (v1[1][j] - v1[0][j]) / dx;
			ksi_2[0][j] = ksi_1[0][j];

			ksi_1[n-1][j] = (u1[n-1][j+1] - u1[n-1][j-1]) / (2*dy) - (v1[n-1][j] - v1[n-2][j]) / dx;
			ksi_2[n-1][j] = ksi_1[n-1][j];
	}
	// psi, ksi ����� � ������
	// 1� �������
	for (i = 0; i < n - 1; i++) {
		ksi[0][i][m-1] = (u1[i][m-1] - u1[i][m-2]) / dy - (v1[i+1][m-1] - v1[i][m-1]) / dx;
		ksi[0][i][0] = (u1[i][1] - u1[i][0]) / dy - (v1[i+1][0] - v1[i][0]) / dx;
	}
	ksi[0][n-1][m-1] = ksi[0][0][m-1]; ksi[0][n-1][0] = ksi[0][0][0]; // ����� ������� �����

	for (i = 1; i < n - 1  ; i++) {
		psi[0][i][m-1] = psi[0][i][m-2] + u1[i][m-1] * dy;
		psi[0][i][0] = psi[0][i][1] - u1[i][0] * dy; 	//�����. ��� � �� ������

		// 2� �������
		//ksi[0][i][m-1] = 2 * (psi[0][i][m-2] - psi[0][i][m-1] + u[0][i][m-1] * dy) / dy / dy;
		//ksi[0][i][0] = 2 * (psi[0][i][1] - psi[0][i][0] - u[0][i][0] * dy) / dy / dy;
	}


	// output.open("FILE.txt");
	// //for (i = 0; i < n; i++){	output << i * dx << "  " << v[0][i][Gran] << endl; } // ����������
	// for (i = 0; i < n; i++){
	// 	for (j = 0; j < m; j++){
	// 			output << x[i] << "  " << y[j] << "  " << ksi[0][i][j]<<endl;
	// 			//output << x[i] << "  " << y[j] <<"  "<< sqrt(u[0][i][j] * u[0][i][j] + v[0][i][j] * v[0][i][j]) << endl;
	// 				//output << x[i] << "  " << y[j] << "  " << u[0][i][j] << " " << v[0][i][j]
	// 				//	<< " "<< sqrt(u[0][i][j] * u[0][i][j] + v[0][i][j] * v[0][i][j]) << endl;
	// 	}
	// }
	// output.close();
	// //fprintf(gnuplotPipe, "plot 'FILE.txt' with lines \n");
	// fprintf(gnuplotPipe, "plot 'FILE.txt' with image \n");
	// //fprintf(gnuplotPipe, "plot '%s' u 1:2:($3*0.1):($4*0.1):5 with vectors head size 0.,0,0 filled lc palette\n", "FILE.txt");
  // fflush(gnuplotPipe);
	// std::cout << "In the begining:" << '\n';
	// std::cout << "Press Enter" << '\n';
	// getchar();


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double L = 1000, dt = 0.001, t0 = L/u0; // ��������� �������
	double alfa, beta, nu, Re, Pe, Ri, Ra, xi, g;
	beta = 0.0001; xi = 0.02; nu = 0.02; g = 10;
	Re = u0 * L / nu; Pe = u0 * L / xi; Ra = g * beta * Th0 * L*L*L / (nu*xi); Ri = Ra / (Re * Pe);

	//f(1000, 100, 1000, Re, Pe, Ri);
	//	f(10000, 10000, 10, Re, Pe, Ri);
	f(50000, 50000, 5, Re, Pe, Ri);

	int s = 0, s_max = 1000000,K, s_image = 5;
	double eps = 1, max;

	while (s < s_max) {
		// ����� ����� ksi
		for (i = 1; i < n - 1; i++) {
			for (j = 1; j < m - 1; j++) {
					if (u[s%2][i][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][i][j] >= 0) { beta = 0; } else { beta = 1; }

					ksi[(s+1)%2][i][j] = ksi[s%2][i][j]
						- dt/dx * u[s%2][i][j] * ((1 - alfa)*(ksi[s%2][i][j] - ksi[s%2][i-1][j]) + alfa*(ksi[s%2][i+1][j] - ksi[s%2][i][j]))
						- dt/dy * v[s%2][i][j] * ((1 - beta)*(ksi[s%2][i][j] - ksi[s%2][i][j-1]) + beta*(ksi[s%2][i][j+1] - ksi[s%2][i][j]))
							+ dt/Re*((ksi[s%2][i+1][j] - 2 * ksi[s%2][i][j] + ksi[s%2][i-1][j])/dx/dx
							+ (ksi[s%2][i][j+1] - 2 * ksi[s%2][i][j] + ksi[s%2][i][j-1])/dy/dy)
								- dt * Ri/(2*dx) * (theta[s%2][i+1][j] - theta[s%2][i-1][j]); }}

		/////// ������� ksi ����� � ������ (������������� ��������� �������).
		for (j = 1; j < m - 1; j++) {
				if (u[s%2][0][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][0][j] >= 0) { beta = 0; } else { beta = 1; }

				ksi[(s+1)%2][0][j] = ksi[s%2][0][j]
					- dt/dx * u[s%2][0][j] * ((1 - alfa)*(ksi[s%2][0][j] - ksi[s%2][n-2][j]) + alfa*(ksi[s%2][1][j] - ksi[s%2][0][j]))
					- dt/dy * v[s%2][0][j] * ((1 - beta)*(ksi[s%2][0][j] - ksi[s%2][0][j-1]) + beta*(ksi[s%2][0][j+1] - ksi[s%2][0][j]))
						+ dt/Re*((ksi[s%2][1][j] - 2 * ksi[s%2][0][j] + ksi[s%2][n-2][j])/dx/dx
						+ (ksi[s%2][0][j+1] - 2 * ksi[s%2][0][j] + ksi[s%2][0][j-1])/dy/dy)
							- dt * Ri/(2*dx) * (theta[s%2][1][j] - theta[s%2][n-2][j]);
				ksi[(s+1)%2][n-1][j] = ksi[(s+1)%2][0][j]; // ������ = �����
			}

		/////////////////// �������. psi
		K = 0;
		do {
			max = 0;
			// ��� ������������ �����. ������
			for (i = 1; i < n - 1; i++){
				for (j = 1; j < m - 1; j++){
						psi[(s+1) % 2][i][j] = (dy*dy*(psi[s%2][i+1][j] + psi[(s+1)%2][i-1][j]) + dx * dx*(psi[s%2][i][j+1] + psi[(s+1)%2][i][j-1])
							- ksi[s%2][i][j] * dx*dx*dy*dy)
							/ (2. * (dy * dy + dx * dx));
			}}
			//������� ������ � ����� ��� psi
			for (i = 0; i < n; i++) {
				psi[(s+1)%2][i][m-1] = psi[(s+1) % 2][i][m-2] + u[s%2][i][m-1] * dy; 	// ������
				psi[(s+1)%2][i][0] = psi[(s+1) % 2][i][1] - u[s%2][i][0] * dy; // �����. ��� � �� ������
			 }
			////// ������� psi ����� � ������ (������������� ��������� �������)
			for (j = 1; j < m - 1; j++){
						psi[(s+1) % 2][0][j] = (dy*dy*(psi[s%2][1][j] + psi[(s+1) % 2][n-2][j]) + dx * dx*(psi[s%2][0][j+1] + psi[(s+1)% 2][0][j-1])
							- ksi[s%2][0][j] * dx*dx*dy*dy)
							/ (2. * (dy*dy + dx * dx));

						psi[(s+1) % 2][n-1][j] = psi[(s+1) % 2][0][j]; // ������ = �����
			}
			// �������� ����������
			for (i = 1; i < n-1; i++){
				for (j = 1; j < m-1; j++){
						max += abs(psi[(s+1) % 2][i][j] - psi[s%2][i][j]);
						psi[s % 2][i][j] = psi[(s + 1) % 2][i][j]; // ����������, � �� �� ������� �� ��������
			}}

			K++;
			if(K>30000) {
				//std::cout << "stop me please!" << '\n';
				break;
			}
			//std::cout << max << '\n';
		//	getchar();
		} while (max > eps);

		///////////////////////// �������� ������ //////////////////////////////////////
		for (i = 1; i < n - 1; i++){
			for (j = 1; j < m - 1; j++){
				// ����������� �� 2� �������� ������
				u[(s+1)%2][i][j] = (psi[(s+1)%2][i][j+1] - psi[(s+1)%2][i][j-1])/2./dy;
			  v[(s+1)%2][i][j] = -(psi[(s+1)%2][i+1][j] - psi[(s+1)%2][i-1][j])/2./dx;
			}
		}

		// �������� (������������� ��������� �������)
		for (j = 1; j < m - 1; j++){
			// ����������� �� 2� �������� ������
			u[(s+1)%2][0][j] = (psi[(s+1)%2][0][j+1] - psi[(s+1)%2][0][j-1])/2./dy;
			u[(s+1)%2][n-1][j] = u[(s+1)%2][0][j]; // ������ = �����

			v[(s+1)%2][0][j] = -(psi[(s+1)%2][1][j] - psi[(s+1)%2][n-2][j])/2./dx;
			v[(s+1)%2][n-1][j] = v[(s+1)%2][0][j]; // ������ = �����

			// ���
			// 	u[(s + 1) % 2][n - 1][j] = u[(s + 1) % 2][n - 2][j];
			// 	u[(s + 1) % 2][0][j] = u[(s + 1) % 2][n-1][j];
			// 	v[(s + 1) % 2][n - 1][j] = v[(s + 1) % 2][n - 2][j];
			// 	v[(s + 1) % 2][0][j] = v[(s + 1) % 2][n-1][j];
		}

		////////////////////// ���ר� ����������� theta
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
		//// ������� �� theta ����� � ������ (������������� ��������� �������)
		for (j = 1; j < m-1; j++) {
			if (u[s%2][0][j] >= 0) { alfa = 0; } else { alfa = 1; } if (v[s%2][0][j] >= 0) { beta = 0; } else { beta = 1; }

			theta[(s+1)%2][0][j] = theta[s%2][0][j]
				- dt/dx * u[s%2][0][j] * ((1 - alfa)*(theta[s%2][0][j] - theta[s%2][n-2][j]) + alfa*(theta[s%2][1][j] - theta[s%2][0][j]))
				- dt/dy * v[s%2][0][j] * ((1 - beta)*(theta[s%2][0][j] - theta[s%2][0][j-1]) + beta*(theta[s%2][0][j+1] - theta[s%2][0][j]))
					+(1/Pe)*dt * ((theta[s%2][1][j] - 2 * theta[s%2][0][j] + theta[s%2][n-2][j]) /dx/dx
						+ (theta[s%2][0][j+1] - 2 * theta[s%2][0][j] + theta[s%2][0][j-1]) /dy/dy);
			theta[(s+1)%2][n-1][j] = theta[(s+1)%2][0][j]; // ������ = �����
		}

		//////////////// ������� ksi � theta ������� ������� � ������.
		// 1� �������
		for (i = 0; i < n - 1; i++) {
			ksi[(s+1)%2][i][m-1] = (u[(s+1)%2][i][m-1] - u[(s+1)%2][i][m-2]) / dy - (v[(s+1)%2][i+1][m-1] - v[(s+1)%2][i][m-1]) / dx;
			ksi[(s+1)%2][i][0] = (u[(s+1)%2][i][1] - u[(s+1)%2][i][0]) / dy - (v[(s+1)%2][i+1][0] - v[(s+1)%2][i][0]) / dx;
	  }
		ksi[(s+1)%2][n-1][m-1] = ksi[(s+1)%2][0][m-1]; ksi[(s+1)%2][n-1][0] = ksi[(s+1)%2][0][0]; // ����� ������� �����

		for (i = 0; i < n; i++) {
			// 2� �������
			// ksi[(s+1)%2][i][m-1] = 2 * (psi[(s+1)%2][i][m-2] - psi[(s+1)%2][i][m-1] + u[(s+1)%2][i][m-1] * dy) / (dy*dy);
			// ksi[(s+1)%2][i][0] = 2 * (psi[(s+1)%2][i][1] - psi[(s+1)%2][i][0] - u[(s+1)%2][i][0] * dy) / (dy*dy); // ����� ??????

			theta[(s+1)%2][i][m - 1] = theta[s%2][i][m - 2] ;
			theta[(s+1)%2][i][0] = theta[s%2][i][1];
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
					//fprintf(gnuplotPipe, "%lf, %lf, %lf\n", x[i],y[j], sqrt(pow(u[(s+1)%2][i][j],2) + pow(v[(s+1)%2][i][j], 2))	);

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
