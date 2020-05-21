// code provided by the authors

#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

void  rlambda11(int Nlambda[], int lambda[], double p[], int V[], int dime[],
        double a[], double logitmean[], int Hdim[], int nMCMC[], double sima[],
        int simlambda[], int simnz[], int simNpop[], int nunit[], int par_up[],
        int prior[], double up_lambda[], double gzeta[], double sigmal[],
        double simp[])
{
    double rU;
    int i, j,l,k,q,old,h,iterMCMC,ksize,totdime;
    int N=Nlambda[0], H=Hdim[0],cumdime[H], NMCMC=nMCMC[0],n=nunit[0];
    int perm[n], permN[100000];
    int clsize[n];
    double psamp[n], pcluster_M[n*H],  prob_M[n*H], pcluster_M_new[n*H],
    prob_M_new[n*H], pcluster_M_old[n*H],  prob_M_old[n*H];
    double logpcond_a_old[n*H], logpcond_a_new[n*H],
    logpcond_aSC_new[n*H],logpcond_aSC_old[n*H],
    a_prop[H*n],logita[H*n],logitadif[H*n],logitadif_prop[H*n];;
    double psampN[100000],maxP;
    double g=gzeta[0],lqratio,cdir,logitmean_prop[H];
    int Nmax=100000;

    GetRNGstate();
    cumdime[0]=0;
    totdime=0;
    for (h=1;h<H;h++) cumdime[h]=cumdime[h-1]+dime[h-1];
    for (h=0;h<H;h++) totdime=totdime+dime[h];
    for (l=0;l<n;l++){
        for(h=0;h<H;h++){
            logita[h*n+l]    = log(a[h*n+l]  / (1-a[h*n+l]));
            logitadif[h*n+l] = logita[h*n+l] - logitmean[h];
        }
    }

    double p_prop[totdime], s_prop;

    /**********************************************************/
    /* CALCOLO DELLE NUMEROSITA' DEI CLUSTER E DEL LORO NUMERO*/
    /**********************************************************/
    for (l=0;l<n;l++) clsize[l]=0;
    for (l=0;l<n;l++) clsize[lambda[l]]=clsize[lambda[l]]+1;
    ksize=0;
    for (l=0;l<n;l++) {if (clsize[l]>0) ksize=ksize+1;}
    /*********************************************************/



    for (iterMCMC=0;iterMCMC<NMCMC;iterMCMC++){

        if (par_up[0]==1){
            /*******************************************/
            /*********** LAMBDA UPDATING ***************/
            /*******************************************/

            if (prior[0]==0){
                for (l=0;l<(H*n);l++) {
                    pcluster_M[l]=1.0;
                    prob_M[l]=1.0;
                }

                /*****************************************************************/
                /* CALCOLO LE PROBABILITA' DEI CLUSTER OSSERVATI */
                /***************************************************************/

                for(l=0;l<n;l++){
                    old=0;
                    for (i=0;i<l;i++) {if (lambda[l]==lambda[i]) {old=1;break; }}

                    if (old==1) {
                        for (h=0;h<H;h++){
                            prob_M[lambda[l]+n*h]=1.0;
                            for ( i=0;i<l;i++){
                                if (lambda[i]==lambda[l]) 
                                    prob_M[lambda[l]+n*h]=prob_M[lambda[l]+n*h]*((V[i+n*h]==V[l+n*h])*(1-a[h*n+lambda[l]])+a[h*n+lambda[l]]*p[cumdime[h]+V[i+n*h] ]);}

                            prob_M[lambda[l]+n*h]     =a[h*n+lambda[l]]*p[cumdime[h]+V[l+n*h]]*pcluster_M[lambda[l]+n*h]+(1-     a[h*n+lambda[l]])*p[V[l+n*h]+cumdime[h]]*prob_M[lambda[l]+n*h]; 
                        } 
                    }
                    else { 
                        for (h=0;h<H;h++) prob_M[lambda[l]+n*h]=p[cumdime[h]+V[l+n*h]];
                    }


                    for (h=0;h<H;h++)  pcluster_M[lambda[l]+n*h]=prob_M[lambda[l]+n*h];     
                }


            }


            for(j=0;j<n;j++){
                if (unif_rand()<up_lambda[0]){

                    /*****************************************************************/
                    /* CALCOLO LA PROBABILITA' DEI CLUSTER OSSERVATI SENZA L'UNITA J*/
                    /* DEVO RICALCOLARE SOLO IL CLUSTER DI LAMBDA[J]*/
                    /***************************************************************/



                    if (prior[0]==0){
                        if (clsize[lambda[j]]>1) {
                            for (h=0;h<H;h++){
                                prob_M[lambda[j]+n*h]=1.0;
                                for ( i=0;i<n;i++){
                                    if (i!=j  && lambda[i]==lambda[j]) 
                                        prob_M[lambda[j]+n*h]=
                                            prob_M[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a[h*n+lambda[j]])+a[h*n+lambda[j]]*p[cumdime[h]+V[i+n*h] ]); 
                                }
                                pcluster_M[lambda[j]+n*h]=
                                    (pcluster_M[lambda[j]+n*h]-(1-a[h*n+lambda[j]])*p[cumdime[h]+V[j+n*h] ]* prob_M[lambda[j]+n*h]) /
                                        (a[h*n+lambda[j]]*p[cumdime[h]+V[j+n*h] ]);
                            }
                        }
                        if (clsize[lambda[j]]==1) {
                            for (h=0;h<H;h++) 
                                pcluster_M[lambda[j]+n*h]=1.0;
                        }
                    }
                    if (clsize[lambda[j]]==1) ksize=ksize-1;
                    for(q=0; q<n; q++) {
                        psamp[q]=1.0;
                        if (prior[0]==0) {
                            for (h=0;h<H;h++){
                                prob_M[q+n*h]=1.0;
                                for ( i=0;i<n;i++){
                                    if (i!=j  && lambda[i]==q)  prob_M[q+n*h]=prob_M[q+n*h]*((V[i+n*h]==V[j+n*h])*(1-a[h*n+q])+a[h*n+q]*p[cumdime[h]+V[i+n*h] ]);
                                }
                                prob_M[q+n*h]=a[h*n+q]*p[cumdime[h]+V[j+n*h]]*pcluster_M[q+n*h]+(1-a[h*n+q])*p[V[j+n*h]+cumdime[h]]*prob_M[q+n*h];
                                psamp[q]=psamp[q]*prob_M[q+n*h]/pcluster_M[q+n*h];
                            }
                            if  (clsize[q]==0 || (clsize[q]==1 && q==lambda[j])) psamp[q]=psamp[q]*(N-ksize+0.0)/(n-ksize+0.0);}
                        if (prior[0]==1){
                            if (clsize[q]==0 ||(clsize[q]==1 && q==lambda[j]))
                                psamp[q]=psamp[q]*(N-ksize+0.0)/(n-ksize+0.0);
                        }
                    }
                    clsize[lambda[j]]=clsize[lambda[j]]-1;
                    /* record element identities */
                    for (i = 0; i < n; i++){ perm[i] = i;}
                    /* sort the probabilities into descending order */
                    revsort(psamp, perm, n);
                    /* compute cumulative probabilities */
                    for (i = 1 ; i < n; i++)
                        psamp[i] += psamp[i - 1];
                    /* sum to 1 the  cumulative probabilities */
                    for (i = 0 ; i < n; i++){
                        psamp[i] = psamp[i]/psamp[n-1];
                    }
                    /* compute the sample */
                    rU = unif_rand();
                    for (l = 0; l < (n-1); l++) {
                        if (rU <= psamp[l])
                            break;
                    }
                    lambda[j]= perm[l];
                    for (h=0;h<H;h++)  pcluster_M[lambda[j]+n*h]=prob_M[lambda[j]+n*h];
                    clsize[lambda[j]]=clsize[lambda[j]]+1;
                    if (clsize[lambda[j]]==1) ksize=ksize+1;
                }
            }
        }

        /*******************************************/
        /***********  UPDATING  THE VECTOR a *******/
        /*******************************************/
        if (par_up[1]==1) {
            for (l=0;l<n;l++){
                for (h=0;h<H;h++) {  
                    logitadif_prop[h*n+l]=rnorm(logitadif[h*n+l],sqrt(sigmal[0])/(fmax(1.0,clsize[l]+0.0)+0.0));
                    a_prop[h*n+l]=exp(logitadif_prop[h*n+l]+logitmean[h])/(1+exp(logitadif_prop[h*n+l]+logitmean[h]));
                    logpcond_a_old[h*n+l]=0.0;
                    logpcond_a_new[h*n+l]=0.0;
                } 
            }
            for(h=0;h<H;h++) {
                pcluster_M_new[h*n+lambda[0]]=p[cumdime[h]+V[n*h]];
                pcluster_M_old[h*n+lambda[0]]=p[cumdime[h]+V[n*h]];
                logpcond_a_new[h*n+lambda[0]]=log(pcluster_M_new[h*n+lambda[0]]);
                logpcond_a_old[h*n+lambda[0]]=log(pcluster_M_old[h*n+lambda[0]]); 
            } 
            for(j=1;j<n;j++){
                old=0;
                for (i=0;i<j;i++) {
                    if (lambda[j]==lambda[i]) {old=1;break; }
                }
                if (old==1){
                    for (h=0;h<H;h++){
                        prob_M_new[lambda[j]+n*h]=1.0;
                        prob_M_old[lambda[j]+n*h]=1.0;
                        for ( i=0;i<j;i++){
                            if (lambda[i]==lambda[j]) {
                                prob_M_new[lambda[j]+n*h]=
                                    prob_M_new[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a_prop[h*n+lambda[j]])+a_prop[h*n+lambda[j]]*p[cumdime[h]+V[i+n*h] ]);
                                prob_M_old[lambda[j]+n*h]=
                                    prob_M_old[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a[h*n+lambda[j]])+a[h*n+lambda[j]]*p[cumdime[h]+V[i+n*h] ]);
                            }  
                        }
                        prob_M_new[lambda[j]+n*h] =
                            a_prop[h*n+lambda[j]]*p[cumdime[h]+V[j+n*h]]*pcluster_M_new[lambda[j]+n*h]+
                             (1-a_prop[h*n+lambda[j]])*p[V[j+n*h]+cumdime[h]]*prob_M_new[lambda[j]+n*h];
                        prob_M_old[lambda[j]+n*h] =
                            a[h*n+lambda[j]]*p[cumdime[h]+V[j+n*h]]*pcluster_M_old[lambda[j]+n*h]+
                            (1-a[h*n+lambda[j]])*p[V[j+n*h]+cumdime[h]]*prob_M_old[lambda[j]+n*h];
                        logpcond_aSC_new[lambda[j]+n*h]=log(prob_M_new[lambda[j]+n*h]/pcluster_M_new[lambda[j]+n*h]);
                        logpcond_aSC_old[lambda[j]+n*h]=log(prob_M_old[lambda[j]+n*h]/pcluster_M_old[lambda[j]+n*h]); 
                    }
                }
                else {
                    for (h=0;h<H;h++){
                        logpcond_aSC_new[lambda[j]+n*h]=log(p[cumdime[h]+V[j+n*h]]);
                        logpcond_aSC_old[lambda[j]+n*h]=log(p[cumdime[h]+V[j+n*h]]);
                        prob_M_new[lambda[j]+n*h]=p[cumdime[h]+V[j+n*h]];
                        prob_M_old[lambda[j]+n*h]=p[cumdime[h]+V[j+n*h]];
                    } 
                }
                for (h=0;h<H;h++)  {
                    pcluster_M_new[lambda[j]+n*h]=prob_M_new[lambda[j]+n*h];
                    pcluster_M_old[lambda[j]+n*h]=prob_M_old[lambda[j]+n*h];
                    logpcond_a_new[h*n+lambda[j]]=logpcond_a_new[h*n+lambda[j]]+logpcond_aSC_new[lambda[j]+n*h];
                    logpcond_a_old[h*n+lambda[j]]=logpcond_a_old[h*n+lambda[j]]+logpcond_aSC_old[lambda[j]+n*h];
                }
            }
            for (l=0;l<n;l++){
                for (h=0;h<H;h++){
                    if (unif_rand()<exp(logpcond_a_new[h*n+l]-logpcond_a_old[h*n+l]+
                                dnorm(logitadif_prop[h*n+l],0,sqrt(sigmal[0]),1)-
                                dnorm(logitadif[h*n+l],     0,sqrt(sigmal[0]),1) )){
                        a[h*n+l]         = a_prop[h*n+l];
                        logitadif[h*n+l] = logitadif_prop[h*n+l];
                    }
                }
            }
            /* aggiornamento di logitmean[h] che ora compare in tutti i cluster osservati*/    
            for (h=0;h<H;h++){
                logitmean_prop[h]=rnorm(logitmean[h],sigmal[1]);
                logpcond_a_old[h]=0;
                logpcond_a_new[h]=0;
            }
            for (l=0;l<n;l++){
                for(h=0;h<H;h++) a_prop[h*n+l]=exp(logitadif[h*n+l]+logitmean_prop[h])/(1+exp(logitadif[h*n+l]+logitmean_prop[h]));
            }

            for(h=0;h<H;h++) {
                pcluster_M_new[h*n+lambda[0]]=p_prop[cumdime[h]+V[n*h]];
                pcluster_M_old[h*n+lambda[0]]=p[cumdime[h]+V[n*h]];
                logpcond_a_new[h]=log(pcluster_M_new[h*n+lambda[0]]);
                logpcond_a_old[h]=log(pcluster_M_old[h*n+lambda[0]]); 
            }
            for(j=1;j<n;j++){
                old=0;
                for (i=0;i<j;i++) {if (lambda[j]==lambda[i]) {old=1;break; }}
                if (old==1){

                    for (h=0;h<H;h++){
                        prob_M_new[lambda[j]+n*h]=1.0;
                        prob_M_old[lambda[j]+n*h]=1.0;
                        for ( i=0;i<j;i++){
                            if (lambda[i]==lambda[j]) {
                                prob_M_new[lambda[j]+n*h]=
                                    prob_M_new[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a_prop[h*n+lambda[j]])+a_prop[h*n+lambda[j]]*p[cumdime[h]+V[i+n*h] ]);
                                prob_M_old[lambda[j]+n*h]=
                                    prob_M_old[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a[h*n+lambda[j]])+a[h*n+lambda[j]]*p[cumdime[h]+V[i+n*h] ]);
                            }  
                        }
                        prob_M_new[lambda[j]+n*h]=
                            a_prop[h*n+lambda[j]]*p[cumdime[h]+V[j+n*h]]*pcluster_M_new[lambda[j]+n*h]+(1-a_prop[h*n+lambda[j]])*p[V[j+n*h]+cumdime[h]]*prob_M_new[lambda[j]+n*h];
                        prob_M_old[lambda[j]+n*h]=
                            a[h*n+lambda[j]]*p[cumdime[h]+V[j+n*h]]*pcluster_M_old[lambda[j]+n*h]+(1-a[h*n+lambda[j]])*p[V[j+n*h]+cumdime[h]]*prob_M_old[lambda[j]+n*h];
                        logpcond_aSC_new[lambda[j]+n*h]=log(prob_M_new[lambda[j]+n*h]/pcluster_M_new[lambda[j]+n*h]);
                        logpcond_aSC_old[lambda[j]+n*h]=log(prob_M_old[lambda[j]+n*h]/pcluster_M_old[lambda[j]+n*h]); 
                    }
                }
                else {
                    for (h=0;h<H;h++){
                        logpcond_aSC_new[lambda[j]+n*h]=log(p[cumdime[h]+V[j+n*h]]);
                        logpcond_aSC_old[lambda[j]+n*h]=log(p[cumdime[h]+V[j+n*h]]);
                        prob_M_new[lambda[j]+n*h]=p[cumdime[h]+V[j+n*h]];
                        prob_M_old[lambda[j]+n*h]=p[cumdime[h]+V[j+n*h]];
                    } 
                }
                for (h=0;h<H;h++)  {
                    pcluster_M_new[lambda[j]+n*h]=prob_M_new[lambda[j]+n*h];
                    pcluster_M_old[lambda[j]+n*h]=prob_M_old[lambda[j]+n*h];
                    logpcond_a_new[h]=logpcond_a_new[h]+logpcond_aSC_new[lambda[j]+n*h];
                    logpcond_a_old[h]=logpcond_a_old[h]+logpcond_aSC_old[lambda[j]+n*h];
                }
            }
            for (h=0;h<H;h++){
                if (unif_rand()<exp(logpcond_a_new[h]-logpcond_a_old[h]+
                            dnorm(logitmean_prop[h],log(sigmal[2]/(1-sigmal[2])),sqrt(sigmal[1]),1)-
                            dnorm(logitmean[h],log(sigmal[2]/(1-sigmal[2])),sqrt(sigmal[1]),1)))
                    logitmean[h]=logitmean_prop[h];}
            for (l=0;l<n;l++){
                for(h=0;h<H;h++) a[h*n+l]=exp(logitadif[h*n+l]+logitmean[h])/(1+exp(logitadif[h*n+l]+logitmean[h]));
            }
        }
        /*******************************************/
        /***********  UPDATING  theta        *******/
        /*******************************************/
        if (par_up[2]==1) {
            cdir=10000.0;
            for (h=0;h<H;h++){
                s_prop=0.0;
                for (k=0;k<dime[h];k++){
                    p_prop[cumdime[h]+k]=rgamma((cdir*p[cumdime[h]+k]+1.0),1.0);
                    s_prop=s_prop+p_prop[cumdime[h]+k];
                }
                for (k=0;k<dime[h];k++) {
                    p_prop[cumdime[h]+k]=p_prop[cumdime[h]+k]/s_prop;
                    /*	 Rprintf("%f %f \n",p_prop[cumdime[h]+k],p[cumdime[h]+k]);*/
                }
                logpcond_a_old[h]=0;
                logpcond_a_new[h]=0;
            }
            for(h=0;h<H;h++) {
                pcluster_M_new[h*n+lambda[0]]=p_prop[cumdime[h]+V[n*h]];
                pcluster_M_old[h*n+lambda[0]]=p[cumdime[h]+V[n*h]];
                logpcond_a_new[h]=log(pcluster_M_new[h*n+lambda[0]]);
                logpcond_a_old[h]=log(pcluster_M_old[h*n+lambda[0]]); 
            }
            for(j=1;j<n;j++){
                old=0;
                for (i=0;i<j;i++) {if (lambda[j]==lambda[i]) {old=1;break; }}
                if (old==1){
                    for (h=0;h<H;h++){
                        prob_M_new[lambda[j]+n*h]=1.0;
                        prob_M_old[lambda[j]+n*h]=1.0;
                        for ( i=0;i<j;i++){
                            if (lambda[i]==lambda[j]) {
                                prob_M_new[lambda[j]+n*h]=
                                    prob_M_new[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a[h*n+lambda[j]])+a[h*n+lambda[j]]*p_prop[cumdime[h]+V[i+n*h] ]);
                                prob_M_old[lambda[j]+n*h]=
                                    prob_M_old[lambda[j]+n*h]*( (V[i+n*h]==V[j+n*h])*(1-a[h*n+lambda[j]])+a[h*n+lambda[j]]*p[cumdime[h]+V[i+n*h] ]);
                            }  
                        }
                        prob_M_new[lambda[j]+n*h]=
                            a[h*n+lambda[j]]*p_prop[cumdime[h]+V[j+n*h]]*pcluster_M_new[lambda[j]+n*h]+
                                (1-a[h*n+lambda[j]])*p_prop[V[j+n*h]+cumdime[h]]*prob_M_new[lambda[j]+n*h];
                        prob_M_old[lambda[j]+n*h]     =
                            a[h*n+lambda[j]]*p[cumdime[h]+V[j+n*h]]*pcluster_M_old[lambda[j]+n*h]+
                                (1-a[h*n+lambda[j]])*p[V[j+n*h]+cumdime[h]]*prob_M_old[lambda[j]+n*h];
                        logpcond_aSC_new[lambda[j]+n*h]=
                            log(prob_M_new[lambda[j]+n*h]/pcluster_M_new[lambda[j]+n*h]);
                        logpcond_aSC_old[lambda[j]+n*h]=
                            log(prob_M_old[lambda[j]+n*h]/pcluster_M_old[lambda[j]+n*h]); 
                    }
                }
                else {
                    for (h=0;h<H;h++){
                        logpcond_aSC_new[lambda[j]+n*h]=log(p_prop[cumdime[h]+V[j+n*h]]);
                        logpcond_aSC_old[lambda[j]+n*h]=log(p[cumdime[h]+V[j+n*h]]);
                        prob_M_new[lambda[j]+n*h]=p_prop[cumdime[h]+V[j+n*h]];
                        prob_M_old[lambda[j]+n*h]=p[cumdime[h]+V[j+n*h]];
                    } 
                }
                for (h=0;h<H;h++)  {
                    pcluster_M_new[lambda[j]+n*h]=prob_M_new[lambda[j]+n*h];
                    pcluster_M_old[lambda[j]+n*h]=prob_M_old[lambda[j]+n*h];
                    logpcond_a_new[h]=logpcond_a_new[h]+logpcond_aSC_new[lambda[j]+n*h];
                    logpcond_a_old[h]=logpcond_a_old[h]+logpcond_aSC_old[lambda[j]+n*h];
                }
            }
            for (h=0;h<H;h++){
                lqratio=0.0;
                for (k=0;k<dime[h];k++) {
                    lqratio=lqratio+(cdir*p[cumdime[h]+k])*log(p_prop[cumdime[h]+k])-(cdir*p_prop[cumdime[h]+k])*log(p[cumdime[h]+k])+
                        lgamma(cdir*p_prop[cumdime[h]+k]+1.0)-lgamma(cdir*p[cumdime[h]+k]+1.0);
                }
                // WARNING: modified lqratio to -lqratio here
                if (unif_rand()<exp(logpcond_a_new[h]-logpcond_a_old[h]-lqratio)){
                    for (k=0;k<dime[h];k++) p[cumdime[h]+k]=p_prop[cumdime[h]+k];
                }
            }
        }    
        if (par_up[3]==1) {
            maxP=-100000.0;
            for(i=0;i<Nmax;i++){
                psampN[i]=lgamma(ksize+i+1)-lgamma(i+1)-(n+g)*log(ksize+i);
                if (psampN[i]>maxP) maxP=psampN[i];
            }
            for(i=0;i<Nmax;i++) psampN[i]=exp(psampN[i]-maxP);
            for (i = 0; i < Nmax; i++){permN[i] = i ;}
            revsort(psampN, permN, Nmax);
            for (i = 1 ; i < Nmax; i++) psampN[i] += psampN[i - 1];
            for (i = 0 ; i < Nmax; i++){ psampN[i] = psampN[i]/psampN[Nmax-1];}
            rU = unif_rand();
            for (l = 0; l < (Nmax-1); l++) {
                if (rU <= psampN[l]) break;
            }
            N=ksize+permN[l];
        }

        simnz[iterMCMC]=ksize;
        simNpop[iterMCMC]=N;
        Rprintf("%d %d %d",iterMCMC, ksize,N);
        for (h=0;h<H;h++){
            sima[H*iterMCMC+h]=exp(logitmean[h])/(1.0+exp(logitmean[h]));
            for (k=0;k<dime[h];k++) simp[totdime*iterMCMC+cumdime[h]+k]=p[cumdime[h]+k];
        }
        for (i=0;i<n;i++) simlambda[n*iterMCMC+i]=lambda[i];
        Rprintf("\n ");
    }
    PutRNGstate();
}


