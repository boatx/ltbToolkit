#include "ltb.h"
#include "hellaby87.h"
#include <string.h>
#include <ltb_integrate.h>
unsigned int N=2000;
unsigned int Nt=5;
double r_min;
double r_hard_min;
double r_hard_max;
//Set to 1 if you want to make plot of radial proper length
int calculate_grr_sq = 1;

void get_r_in(gsl_vector *rVec, int N)
{
    gsl_vector_gen(rVec,N,r_hard_min+3*r_min,r_hard_max - 3*r_min);
}

void get_t_in(gsl_vector *tVec, gsl_vector *tB, double tC, double t_0)
{
    double min_tB = gsl_vector_min(tB);
    double max_tB = gsl_vector_max(tB);
    tVec->data[0]=min_tB+0.5*(max_tB-min_tB);
    tVec->data[1]=min_tB+1.05*(max_tB-min_tB);
    tVec->data[2]=0.5*t_0;
    tVec->data[3]=1.0*t_0;
    tVec->data[4]=tC-2.0*(max_tB-min_tB);
}

double get_t0(double R_C)
{
    //tB(\pi) = 0
    //calculate t0 = 0.5*H_0^2R_C^3*xi_0

    double cos_eta, eta_0, xi_0;
    cos_eta = 1.0 - 2.0/(R_C*R_C*omm*H_0*H_0);
    eta_0 = acos(cos_eta);
    xi_0 = eta_0 - sin(eta_0);
    return 0.5*omm*H_0*H_0*R_C*R_C*R_C*xi_0;
}


int main(void)
{
    tB_norm=1;

    unsigned int i=0;
    Ltb_model *model = ltb_model_alloc(N,Nt,&M_hellaby87,&E_hellaby87,&tB_hellaby87);

    //curvature of universe
    double R_C = 1.0/(H_0*sqrt(omm-1.0));
    double t_0  = get_t0(R_C);

    //initialize global variables from
    r_min=1e-6;
    if (calculate_grr_sq)
    {
        r_hard_min = 0.0;
        r_hard_max = pi;
    }
    else
    {
        r_hard_min=-3.1*pi;
        r_hard_max= 3.1*pi;
    }


    printf("r_hard_min=%g\n",r_hard_min);
    printf("r_hard_max=%g\n",r_hard_max);
    printf("r_min=%g\n",r_min);

    //generate r
    get_r_in(model->r_vec,N);

    //set values of function M,E,tB
    ltb_model_set_fun(model);

    //generate tVec
    //time of big crunch
    double tC = tC_hellaby87(model->M_vec,model->E_vec);
    get_t_in(model->t_vec, model->tB_vec,tC,t_0);
    gsl_write_two_vec(model->t_vec,model->t_vec,"t.dat");
    ltb_model_set_fun(model);
    printf("cur=%d\n",model->cur);
    printf("tC=%g\n",tC);

    //allocate Interp_data for spline interpolation of \eta(\xi)
    Interp_data *eta_interp=interp_data_alloc_for_ltb(model,eta_n);
    ltb_model_write_fun(model,"model_fun.dat");
    model->R_mat = gsl_matrix_alloc(Nt,N);
    //calculate areal radius
    ltb_model_set_R(model,eta_interp);
    gsl_write_mat(model->R_mat,"R.dat",NULL);

    //first derivatives

    //M_r analytic
    puts("M_r");
    gsl_vector *M_r_vec = gsl_vector_alloc(model->r_size);
    for(i=0;i<M_r_vec->size;++i)
    {
        M_r_vec->data[i]=M_r_hellaby87(model->r_vec->data[i]);
    }
    gsl_write_two_vec(model->r_vec,M_r_vec,"M_r.dat");

    //E_r analytic
    puts("E_r");
    gsl_vector *E_r_vec = gsl_vector_alloc(model->r_size);
    for(i=0;i<E_r_vec->size;++i)
    {
        E_r_vec->data[i]=E_r_hellaby87(model->r_vec->data[i]);
    }
    gsl_write_two_vec(model->r_vec,E_r_vec,"E_r_a.dat");

    //tB_r analytic
    gsl_vector *tB_r_vec = gsl_vector_alloc(model->r_size);
    for (i=0;i<model->r_size;++i)
    {
        tB_r_vec->data[i] = tB_r_hellaby87(model->r_vec->data[i]);
    }
    gsl_write_two_vec(model->r_vec,tB_r_vec,"tB_r.dat");

    //R_t
    puts("R_t");
    gsl_matrix *R_t = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_R_t(model,eta_interp,R_t);
    gsl_write_mat(R_t,"R_t.dat",NULL);

    //R_r
    puts("R_r");
    gsl_matrix *R_r = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_R_r(model,R_t,E_r_vec,M_r_vec,tB_r_vec, R_r);
    gsl_write_mat(R_r,"R_r.dat",NULL);

    //R_rt
    puts("R_rt");
    gsl_matrix *R_rt = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_R_rt(model,R_r,R_t,E_r_vec,M_r_vec,R_rt);
    gsl_write_mat(R_rt,"R_rt.dat",NULL);

    //rho
    puts("rho");
    gsl_matrix *rho = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_rho(model,M_r_vec,R_r,rho);
    gsl_write_mat(rho,"rho.dat",NULL);

    //H
    puts("Hubble parameter");
    gsl_matrix *exp = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_exp(model, R_r, R_t, R_rt,exp);
    gsl_write_mat(exp,"exp.dat",NULL);

    //Ricci
    puts("Ricci tensor");
    gsl_matrix *ric = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_Ricci(model, E_r_vec,R_r,ric);
    gsl_write_mat(ric,"ric.dat",NULL);

    if (calculate_grr_sq)
    {
        puts("grr");
        gsl_matrix *grr = gsl_matrix_alloc(model->t_size,model->r_size);
        ltb_model_grr_sq(model->E_vec,R_r,grr);
        gsl_write_mat(grr,"grr.dat",NULL);

        gsl_matrix *integ = gsl_matrix_alloc(model->t_size,model->r_size);

        double dr = model->r_vec->data[1] - model->r_vec->data[0];
        //integrate grr_sq to get radial proper length
        //1 - ommit first column
        gsl_matrix_integrate_trap_range(grr,
                                        dr,
                                        1,model->r_size,
                                        0,model->t_size,
                                        integ);


        gsl_write_mat(integ,"integr.dat",NULL);
        gsl_matrix_free(integ);
        gsl_matrix_free(grr);
    }
    //freeing
    gsl_matrix_free(ric);
    gsl_matrix_free(exp);
    gsl_matrix_free(rho);
    gsl_matrix_free(R_r);
    gsl_matrix_free(R_t);
    gsl_matrix_free(R_rt);
    gsl_vector_free(E_r_vec);
    gsl_vector_free(M_r_vec);
    gsl_vector_free(tB_r_vec);
    interp_data_free(eta_interp);
    ltb_model_free(model);
    return 0;
}
