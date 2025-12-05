typedef	struct {double norme;
						double num;
						double den;
	} residu_stats;
residu_stats norme_residu(int n,int *ia, int *ja, double *a, double beta2, double *phi);
