#include "flrw.h"
#include <ltb.h>
#include <string.h>
#include <omp.h>

int main(int argc, char** argv)
{

    unsigned int Nt=1000;
    unsigned int  N=1000;
    unsigned int i=0;
    unsigned int nt=4;
    omp_set_num_threads(nt);
    chunk_size = Nt*N/nt;
    tB_norm=0;
    Ltb_model *model = ltb_model_alloc(N,Nt,&M_flrw,&E_flrw,&tB_flrw);

    if (argc >= 2)
    {
        printf("a=%s\n",argv[1]);
        initialise_flrw((double)atof(argv[1]),0.0,0.0);
    }
    else
    {
        //initialise parameter for flrw model
        //initialise_flrw(0.5,0.0,0.0);
        initialise_flrw(1.5,0.0,0.0);
        //initialise_flrw(1.0,0.0,0.0);
    }

    printf("number of threads=%d\n",nt);

    double r_min=1e-6;
    double r_hard_min= 0.1;
    double r_hard_max= 0.8;

    printf("r_hard_min=%g\n",r_hard_min);
    printf("r_hard_max=%g\n",r_hard_max);
    printf("r_min=%g\n",r_min);

    //generate r
    gsl_vector_gen(model->r_vec,N,r_hard_max + 3.0*r_min,r_hard_max-3.0*r_min);
    ltb_model_set_fun(model);

    //generate tVec
    //Om_m=1.0;
    //double t_age=6.518630472/Gyr_per_Gpc;
    //Om_m=1.5;
    double t_age=5.976859613/Gyr_per_Gpc;
    //Om_m=0.5;
    //double t_age=7.368166820/Gyr_per_Gpc;

    double t_min=1e-5;
    double t_max=t_age;
    printf("t_min=%g t_max=%g\n",t_min,t_max);

    print_param_flrw(t_min,t_max);
    //export parameters for python
    export_param_flrw("flrw_param.py",t_min,t_max);
    export_phys_const("phys_const.py");

    gsl_vector_gen(model->t_vec,Nt,t_min,t_max);
    gsl_write_two_vec(model->t_vec,model->t_vec,"t.dat");

    printf("cur=%d\n",model->cur);
    Interp_data *eta_interp = interp_data_alloc_for_ltb(model, eta_n);

    ltb_model_write_fun(model,"model_fun.dat");
    model->R_mat = gsl_matrix_alloc(Nt,N);
    ltb_model_set_R(model,eta_interp);
    //set_R(model->t_vec, model->tB_vec, model->M_vec, model->E_vec, model->R_mat, model->cur);
    gsl_write_mat(model->R_mat,"R.dat",NULL);


    //M_r
    puts("Calculating M_r");
    gsl_vector *M_r_a_vec = gsl_vector_alloc(model->r_size);
    for(i=0;i<M_r_a_vec->size;++i)
    {
        M_r_a_vec->data[i]=M_r_flrw(model->r_vec->data[i]);
    }

    //E_r
    puts("Calculating E_r");
    gsl_vector *E_r_a_vec = gsl_vector_alloc(model->r_size);
    for(i=0;i<E_r_a_vec->size;++i)
    {
        E_r_a_vec->data[i]=E_r_flrw(model->r_vec->data[i]);
    }

    //tB_r
    puts("Calculating tB_r");
    gsl_vector *tB_r_a = gsl_vector_alloc(model->r_size);
    for (i=0;i<model->r_size;++i)
    {
        tB_r_a->data[i] = tB_r_flrw(model->r_vec->data[i]);
    }

    gsl_matrix *R_t_a = gsl_matrix_alloc(model->t_size,model->r_size);
    gsl_matrix *R_rt_a = gsl_matrix_alloc(model->t_size,model->r_size);
    gsl_matrix *R_r_a = gsl_matrix_alloc(model->t_size,model->r_size);

    //R_t
    puts("R_t");
    ltb_model_R_t(model,eta_interp,R_t_a);
    gsl_write_mat(R_t_a,"R_t_a.dat",NULL);

    //R_r
    puts("R_r");
    ltb_model_R_r(model,R_t_a,E_r_a_vec,M_r_a_vec,tB_r_a, R_r_a);
    gsl_write_mat(R_r_a,"R_r_a.dat",NULL);

    //R_rt
    puts("R_rt");
    ltb_model_R_rt(model,R_r_a,R_t_a,E_r_a_vec,M_r_a_vec,R_rt_a);
    gsl_write_mat(R_rt_a,"R_rt_a.dat",NULL);


    //rho
    puts("rho");
    gsl_matrix *rho_mat_a = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_rho(model, M_r_a_vec,R_r_a,rho_mat_a);
    gsl_write_mat(rho_mat_a,"rho_a.dat",NULL);

    //H
    puts("H");
    gsl_matrix *exp_a = gsl_matrix_alloc(model->t_size,model->r_size);
    ltb_model_exp(model, R_r_a, R_t_a, R_rt_a,exp_a);
    gsl_write_mat(exp_a,"exp_a.dat",NULL);

    //freeing
    puts("freeing");
    gsl_matrix_free(exp_a);
    gsl_matrix_free(rho_mat_a);
    gsl_matrix_free(R_rt_a);
    gsl_matrix_free(R_r_a);
    gsl_matrix_free(R_t_a);
    gsl_vector_free(tB_r_a);
    gsl_vector_free(E_r_a_vec);
    gsl_vector_free(M_r_a_vec);
    interp_data_free(eta_interp);
    ltb_model_free(model);
    return 0;
}
