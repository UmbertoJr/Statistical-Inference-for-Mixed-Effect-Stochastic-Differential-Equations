import numpy as np
from scipy.stats import norm, multivariate_normal

class dati:
    def __init__(self):
        data = np.genfromtxt('FreyerSutherland.dat')
        self.names = data[1:]
        data = data[1:]
        self.time = data[:,0]
        self.yobs = data[:,1]
        self.subj = data[:,2]

    def select_subj(self, num):
        self.timei = self.time[self.subj == num]
        self.yobsi = self.yobs[self.subj == num]

class Sequential_Monte_Carlo(dati):
    def __init__(self,  number_of_particles = 25, proposal = []):  #IDS, LDS, OD,
        #self.initial_distr_sampler = IDS     # this is the function used to sample the initial distribution
        #self.latent_transition_sampler = LDS   # this samples from transition distribution is used as proposal
        #self.observed_distr = OD   # the observed error is the function used to reweight the particles
        dati.__init__(self)
        self.np = number_of_particles
        if proposal:   # in caso ci siano proposals
            self.init_prop = proposal[0]
            self.transition_prop = proposal[1]
    
    def fit_SMC(self, soggetto , theta, phi=[]):
        self.soggetto = soggetto
        self.select_subj(self.soggetto)
        self.Y = self.yobsi
        self.marginal_particle = [0 for _ in range(len(self.Y))]
        self.marginal_y = 0
        self.theta = theta
        if len(phi)==0:
            print("non ci sono le phii")
            self.phi = theta[:2]
        else:
            self.phi = phi
        self.X_tot = np.zeros((self.np, len(self.Y)))
        self.first_step()  # inizializza x(0)
        self.second_step() # completa le x_i per i = 1,...,len(Y)
    
    def first_step(self):
        X_0 = [self.initial_distr_sampler() for _ in range(self.np)]
        self.X_tot[:,0] = np.array(X_0)
        w_n = [self.observed_distr(x_0, 0) for x_0 in X_0]
        s_w =sum(w_n)
        self.marginal_particle[0] = s_w   # /self.np  ho deciso di toglierlo poiche nel MH scompare
        W_n = np.array(w_n)/s_w
        self.W = W_n
    
    def second_step(self):
        indeces = np.arange(self.np)
        A_n_m = np.random.choice(indeces , size=self.np, p= self.W)
        X_n_m = self.X_tot[A_n_m.tolist(),0]
        self.X_tot[:,0] = X_n_m
        for t in range(1,len(self.Y)):
            X_n = [self.latent_transition_sampler(x_1,t) for x_1 in X_n_m]
            self.X_tot[:,t]= np.array(X_n)
            
            w_n = [self.observed_distr(x, t) for x in X_n]
            s_w =sum(w_n)
            
            self.marginal_particle[t] = s_w    #/ self.np  ho deciso di toglierlo poiche nel MH scompare
            #print(t, self.marginal[t], s_w)
            self.W = np.array(w_n)/s_w
            
            A_n_m = np.random.choice(indeces, size=self.np, p= self.W)
            
            self.X_tot[indeces.tolist(), 0:t] = self.X_tot[A_n_m.tolist(),0 : t]
            X_n_m = self.X_tot[:,t]
            
        self.marginal_y = np.prod(np.array(self.marginal_particle))
        
    def sampling_a_path(self):
        indeces = np.arange(self.np)
        part = np.random.choice(indeces, size=1, p= self.W)
        return self.X_tot[part, :]
    
    def initial_distr_sampler(self):
        return 1e-7
    
    def observed_distr(self, x_, t_ ):
        #foo = (x_ - self.Y[t_])/ self.theta[5]
        #p = 1/(np.sqrt(2*np.pi)* self.theta[5]) * np.exp(-(foo)**2/2)
        p = norm.pdf(x_,self.Y[t_], self.theta[5])
        if p ==0:
            p = 1e-100
        #print(t_, x_, p,p_s,  x_ - self.Y[t_], self.theta[5], foo)
        
        return p
    
    def latent_transition_sampler(self, x_t, time):
        ###   beta = theta[2]
        t_j = self.timei[time]; t_j_1 = self.timei[time-1]
        deltatime = t_j - t_j_1
        u_j = self.EulerStep(x_t, t_j_1, t_j, deltatime)
        s_j = self.theta[4] * u_j * np.sqrt(deltatime)
        x_t1 = u_j + s_j* np.random.normal()
        if x_t1 < 1e-30:
            x_t1 = 1e-30
        
        return x_t1
        
        
    def EulerStep(self, x_t, time0, time1, deltat):
        x_t1 = x_t
        time = time0
        while time < time1:
            if (time1-time)>deltat:
                dt = deltat
            else:
                dt = time1-time
            x_t1 += self.ModelDrift(x_t1)
            time += dt
        return x_t1
    
    def ModelDrift(self, X):
        return self.phi[0]* np.log(self.phi[1]/X)* X
        
        
        
        

class Particle_marginal_Metropolis_Hastings(Sequential_Monte_Carlo):
    def __init__(self, number_of_iterations=100):
        self.NI = number_of_iterations
        #self.SMC = Sequential_Monte_Carlo()
        Sequential_Monte_Carlo.__init__(self)
      
    
    def sample_PMCMC(self, soggetto, theta):
        self.select_subj(soggetto)
        self.X_PMCMC = np.zeros((self.NI, len(self.yobsi)))
        self.phi_PMCMC = np.zeros((self.NI, 2))
        self.marginal_pmcmc = np.zeros(self.NI)
        
        self.theta = theta
        
        #inizializza la prior... faccio questo per velocizzare
        O = np.diag([self.theta[2]**2,self.theta[3]**2])
        self.prior = multivariate_normal(self.theta[:2], O)
        
        accept = 0
        
        phi = self.proposalPhi_sample(np.array(theta[:2]))  # vedere se cambiare...  
        self.phi_PMCMC[0,:] = phi
        
        self.fit_SMC(  soggetto , theta, phi  )
        X = self.sampling_a_path()
        self.X_PMCMC[0, :] = X
        self.marginal_pmcmc[0] = self.marginal_y
        
        
        for i in range(1, self.NI):
            phi_star = self.proposalPhi_sample(phi)
            
            self.fit_SMC( soggetto, theta , phi_star)
            X_star = self.sampling_a_path();
            marginal_star = self.marginal_y
            
            num = (marginal_star* self.prior_phi(phi_star) * self.proposal_evaluation(phi_star, phi))
            den = (self.marginal_pmcmc[i-1]* self.prior_phi(phi) * self.proposal_evaluation(phi, phi_star) )
            prob_acceptance = np.min([1. , num/den ])
            
            if prob_acceptance > np.random.uniform():
                accept += 1
                self.X_PMCMC[i,:] = X_star
                self.phi_PMCMC[i,:] = phi_star
                self.marginal_pmcmc[i] = marginal_star
                X = X_star
                phi = phi_star
            else:
                self.X_PMCMC[i,:] = X
                self.phi_PMCMC[i,:] = phi
                self.marginal_pmcmc[i] = self.marginal_pmcmc[i-1]

            
        print("acceptance rate = ", accept/self.NI)
        return [phi, X]
    
    
    def proposalPhi_sample(self, prec_phi):
        deviance = 0.1*prec_phi
        deviance[deviance < 1e-12] = 1e-12
        sampled = np.array([prec_phi[i] + deviance[i]*np.random.normal() for i in range(2)])
        sampled[sampled < 1e-12] = 1e-12
        return sampled
    
    def prior_phi(self, phi ):
        return self.prior.pdf(phi)
    
    def proposal_evaluation(self, phi_new, phi_mean):
        var = np.array((0.1 * phi_mean)**2)
        var[var < 1e-12] = 1e-12
        variance = np.diag(var)
        #print(phi_new, phi_mean, multivariate_normal( phi_mean, variance).pdf(phi_new))
        return multivariate_normal( phi_mean, variance).pdf(phi_new)