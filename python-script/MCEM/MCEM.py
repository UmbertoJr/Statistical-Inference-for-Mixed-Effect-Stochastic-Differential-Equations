from particle_markov_chain_monte_carlo import Particle_marginal_Metropolis_Hastings, dati
from scipy.optimize import minimize,  differential_evolution, fmin
from scipy.stats import norm, multivariate_normal
import numpy as np

class Monte_Carlo_Expectation_Maximization(dati):
    def __init__(self, turni_da_ottimizzare = 30, campioni_monte_carlo = 10):
        self.M = turni_da_ottimizzare
        self.MC = campioni_monte_carlo
        self.theta_history = np.zeros((6,self.M))
        self.log_lik_history = np.zeros(self.M)
        
        self.PMCMC = Particle_marginal_Metropolis_Hastings(number_of_iterations=300)  #number_of_iterations=1000
        
        dati.__init__(self); 
        
        self.PhiALL = np.zeros((self.MC, 6, 2))
        self.XALL = np.zeros((self.MC, 6, 18))
        
        self.PARMIN = [0.01, 1e-5, 0.002, 1e-4, 0.0001, 0.3e-5] ; self.PARMAX = [1, 1e-1, 0.5, 0.5, 0.1, 5e-5]
        self.bounds = []
        for i in range(6):
            self.bounds.append(tuple([self.PARMIN[i], self.PARMAX[i]]))
        
        
    def run(self, theta_init):
        self.theta_choose = theta_init
        self.phi_subject = np.array(theta_init[:2])
        
        for em in range(self.M):
            print("\n### EM iterazione : ", em)
            for iterazione in range(0, self.MC):
                for subject in range(1,7):
                    if iterazione > 0:
                        self.phi_subject = self.PhiALL[iterazione-1, subject-1, :]
                    phi, X = self.PMCMC.sample_PMCMC(subject, self.theta_choose, self.phi_subject)
                    self.PhiALL[iterazione, subject-1, :] = phi  #self.PMCMC.phi_PMCMC[-1,:] ; 

                    len_subj = len(X)
                    self.XALL[iterazione, subject-1, :len_subj] = X   #self.PMCMC.X_PMCMC[-1,:]



            self.theta_choose = self.ottimizza_MC_approx()
            self.theta_history[:,em] = self.theta_choose
            self.log_lik_history[em] = self.MonteCarlo_Approximation(self.theta_choose)
            
            #if iterazione==0:
             #   self.PMCMC = Particle_marginal_Metropolis_Hastings(number_of_iterations=100)
                
        return self.theta_choose
    
    
    def ottimizza_MC_approx(self):
        
        #                self.iterazione = iterazione
        res = minimize(self.MonteCarlo_Approximation , self.theta_choose, method="Nelder-Mead")
                       
        print(res)
        return res.x
    
    def MonteCarlo_Approximation(self, theta ):
        for i in range(len(theta)):
            if theta[i] > self.PARMAX[i] or theta[i] < self.PARMIN[i]:
                return 1e+300
        MC_log_l = 0
        for mc_ in range(self.MC):
            self.iterazione = mc_
            logPY = self.Compute_LogLikelihood(theta)
            MC_log_l += logPY
        return Qm_1 + self.lambda_ * (logPY - Qm_1)
        
    def Compute_LogLikelihood(self, theta):
        logpY = 0
        #print(theta)
        for i in range(len(theta)):
            if theta[i] > self.PARMAX[i] or theta[i] < self.PARMIN[i]:
                #print(i, theta[i], self.PARMAX[i], self.PARMIN[i])
                return 1e+300
        try:
            self.phii_distr = multivariate_normal(theta[:2], np.diag(np.array(theta[2:4])**2))
        except:
            print("problema")
            print(theta)
        for subj_ in range(1,7):
            self.select_subj(subj_)

            self.phi_i = self.PhiALL[self.iterazione,subj_ - 1, :]
            self.X_i = self.XALL[self.iterazione, subj_ - 1, :len(self.timei)]
            
            pYigivenXi   = self.ModelCondDensityYgivenXevaluate(self.yobsi[0], self.X_i[0], theta);  
            pXigivenPhii = self.ModelTransitionDensityXevaluate_initial(theta)
            for i_ in range(1, len(self.timei)):
                
                pYijgivenXij = self.ModelCondDensityYgivenXevaluate(self.yobsi[i_], self.X_i[i_], theta) 
                pYigivenXi   *=  pYijgivenXij
                #print(pYigivenXi, pYijgivenXij)
                
                pXijgivenPhii = self.ModelTransitionDensityXevaluate(i_, theta)
                pXigivenPhii = pXigivenPhii * pXijgivenPhii;
                
                
            if pYigivenXi==0:
                pYigivenXi=1e-300
            
            if pXigivenPhii==0:
                pXigivenPhii = 1e-300
                
                
            logpYigivenXi = np.log(pYigivenXi)
            logpXigivenPhii = np.log(pXigivenPhii)     
            logpPhii = np.log(self.ModelAprioriDensityPhievaluate(theta))
            
            logpY += (logpYigivenXi + logpXigivenPhii + logpPhii)
            #print(subj_, logpY, logpYigivenXi, logpXigivenPhii, logpPhii)
            #print("\n\n")
        return -logpY
    
    
    def PreviousQm(self, theta):        
        current_iteration = self.iterazione
        self.iterazione = 0
        Qm = self.Compute_LogLikelihood(theta)
        for it_ in range(1, current_iteration):
            self.iterazione = it_
            logPY = self.Compute_LogLikelihood(theta)
            Qm += self.lambda_ * (logPY - Qm)
        return Qm
    
    def ModelAprioriDensityPhievaluate(self, theta):
        q = self.phii_distr.pdf(self.phi_i)
        if q==0:
            q = 1e-300
        return q
    
    def ModelCondDensityYgivenXevaluate(self, y_, x_, theta):
        p = norm.pdf(y_, x_, theta[5])
        if p==0: 
            p = 1e-300
        return p
    
    def ModelTransitionDensityXevaluate_initial(self, theta):
        beta = theta[4]
        time1 = self.timei[0]; time0 = 0
        x_t1 = 1e-7
        time = time0
        deltat =  (time1 - time0)/ 100
        while time < time1:
            if (time1-time)>deltat:
                dt = deltat
            else:
                dt = time1-time
            x_t1 += self.ModelDrift(x_t1 )* dt
            time += dt
        uj = x_t1
        sj = beta * uj * np.sqrt(self.timei[0])

        p = norm.pdf(self.X_i[0] ,  uj, sj)
        if p == 0:
            p = 1e-300      
            
        #print(p, uj, sj)
        return p
    
 
        
    
    def ModelTransitionDensityXevaluate(self, time_, theta):
        beta = theta[4]
        
        uj = self.EulerStep(time_)
        sj = beta * uj * np.sqrt(self.timei[time_] - self.timei[time_- 1])

        p = norm.pdf(self.X_i[time_] ,  uj, sj)
        if p == 0:
            p = 1e-300      
            
        #print(p, uj, sj)
        return p
    
    def EulerStep(self, time_):
        time1 = self.timei[time_]; time0 = self.timei[time_ - 1]
        x_t1 = self.X_i[time_ -1]
        time = time0
        deltat =  (time1 - time0)/ 100
        while time < time1:
            if (time1-time)>deltat:
                dt = deltat
            else:
                dt = time1-time
            x_t1 += self.ModelDrift(x_t1 )* dt
            time += dt
        return x_t1
    
    def ModelDrift(self, X):
        return self.phi_i[0]* np.log(self.phi_i[1]/X)* X
    