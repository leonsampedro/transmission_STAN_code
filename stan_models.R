library(rstan)


#list for stan input (each model requires the specific data)
data.R = list(
  # patients at risk
  atRisk = database$patients_at_risk,
  # binary outcome (acquire Kpn or Eco with pOXA-48)
  events = acquire_pOXA48,
  # patients already colonised by day
  colonised = database$patients_colonised_by_day,
  # number of wards
  N_ward = 4,
  # index of wards
  ward = database$index_of_wards,
  # index of intervals
  interval_index = database$interval_index,
  # number of total intervals (max)
  N_intervals = max(database$interval_index),
  # interval length
  interval_length = interval_length,
  # number of observations
  N=nrow(database),
  # binary previously acquired pOXA-48 by other species
  previously_colonised = database$previously_colonised
)

# MODEL 1

MODEL_1 <- '
data{
int<lower=1> N;                               //number of observations
int<lower=1> N_intervals;                     //total number of intervals
int<lower=0, upper=1> events[N_intervals];    //binary outcome (acquire klebsiella)
int<lower=0> colonised[N];                    //patients already colonised
int<lower=1> N_ward;                          //number of wards (4)
int<lower=1> ward[N];                         //references the ward
int<lower=1> interval_index[N];               //index of intervals
int<lower=1> interval_length[N_intervals];    //interval length
}
parameters{
real alpha;                                   //alpha (intercept), now indexed by ward
}
model{
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
int counter;
alpha ~ normal(0, 10);                        //alpha (intercept)

for(i in 1:N){                                //iterate through patient days (observations)

real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
events ~ bernoulli(interval_prob);              //regression
}
generated quantities{
//calculate log likelihood of model - allows for model comparison
vector[N_intervals] interval_prob;              //vector to store interval probabilities
vector[N] p;                                    //vector to store daily probabilities
int counter;
vector[N_intervals] log_lik;                    //stores log-likelihood of the model

for(i in 1:N){                                  //iterate through patient days (observations)

real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha);
p[i] = day_odds / (day_odds + 1); 
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
//calculate log-likelihood of model
for(k in 1:N_intervals){
log_lik[k] = bernoulli_lpmf(events[k] | interval_prob[k]);
}

}
'


# KPN MODEL
# run the model - note we will now run three parallel chains
RESULT_MODEL_5 <- stan(model_code=MODEL_5, data=data.R,
                              iter=3500,warmup=1500,chains=4, cores=4, thin=2,
                              init_r=8, pars=c("alpha", "beta", "mu", "gamma", "sigma", "log_lik"))


# MODEL 2

MODEL_2 <- '
data{
int<lower=1> N;                               //number of observations
int<lower=1> N_intervals;                     //total number of intervals
int<lower=0, upper=1> events[N_intervals];    //binary outcome (acquire klebsiella)
int<lower=0> colonised[N];                    //patients already colonised
int<lower=1> N_ward;                          //number of wards (4)
int<lower=1> ward[N];                         //references the ward
int<lower=1> interval_index[N];               //index of intervals
int<lower=1> interval_length[N_intervals];    //interval length
}
parameters{
real alpha;                                   //alpha (intercept), now indexed by ward
real beta;                                    //beta (coefficient)
}
model{
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
int counter;
alpha ~ normal(0, 10);                        //alpha (intercept)
beta ~ normal(0, 10);                         //prior for beta (slope, linked to force of infection)

for(i in 1:N){                                //iterate through patient days (observations)

real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta*colonised[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
events ~ bernoulli(interval_prob);              //regression
}
generated quantities{
//calculate log likelihood of model - allows for model comparison
vector[N_intervals] interval_prob;              //vector to store interval probabilities
vector[N] p;                                    //vector to store daily probabilities
int counter;
vector[N_intervals] log_lik;                    //stores log-likelihood of the model

for(i in 1:N){                                  //iterate through patient days (observations)

real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta*colonised[i]);
p[i] = day_odds / (day_odds + 1); 
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
//calculate log-likelihood of model
for(k in 1:N_intervals){
log_lik[k] = bernoulli_lpmf(events[k] | interval_prob[k]);
}

}
'





# MODEL 3
MODEL_3 <-'
data{
int<lower=1> N;                               //number of observations
int<lower=1> N_intervals;                     //total number of intervals
int<lower=0, upper=1> events[N_intervals];    //binary outcome (acquire klebsiella)
int<lower=0> colonised[N];                    //patients already colonised
int<lower=1> N_ward;                          //number of wards (4)
int<lower=1> ward[N];                         //references the ward
int<lower=1> interval_index[N];               //index of intervals
int<lower=1> interval_length[N_intervals];    //interval length
}
parameters{
real alpha;                                   //alpha (intercept), now indexed by ward
real beta[N_ward];                            //beta (coefficient)
real mu;                                      //global mean of intercepts
real<lower=0> sigma;                          //variance of intercepts
}
model{
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
int counter;
alpha ~ normal(0, 10);                        //alpha (intercept)
beta ~ normal(mu, sigma);                     //prior for beta (slope, linked to force of infection)
mu ~ normal(0, 5);                            //prior for global mean
sigma ~ normal(0, 5);                         //prior for variance

for(i in 1:N){                                //iterate through patient days (observations)
real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta[ward[i]]*colonised[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
events ~ bernoulli(interval_prob);              //regression
}
generated quantities{
//calculate log likelihood of model - allows for model comparison
vector[N_intervals] interval_prob;                //vector to store interval probabilities
vector[N] p;                                      //vector to store daily probabilities
vector[N_intervals] log_lik;                      //stores log-likelihood of the model
int counter;

for(i in 1:N){                                    //iterate through patient days (observations)
real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta[ward[i]]*colonised[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
//calculate log-likelihood of model
for(k in 1:N_intervals){
log_lik[k] = bernoulli_lpmf(events[k]|interval_prob[k]);
}
}'






# MODEL 4
MODEL_4 <-'
data{
int<lower=1> N;                                           //number of observations
int<lower=1> N_intervals;                                 //total number of intervals
int<lower=0, upper=1> events[N_intervals];                //binary outcome (acquire klebsiella)
int<lower=0> colonised[N];                                //patients already colonised
int<lower=1> N_ward;                                      //number of wards (4)
int<lower=1> ward[N];                                     //references the ward
int<lower=1> interval_index[N];                           //index of intervals
int<lower=1> interval_length[N_intervals];                //interval length
int<lower=0, upper=1> previously_colonised[N];            //acquire ecoli (binary)
}
parameters{
real alpha;                                           //alpha (intercept), now indexed by ward
real beta;                                            //beta (coefficient)
real gamma;
real mu;                                              //global mean of intercepts
real<lower=0> sigma;                                  //variance of intercepts
}
model{
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
int counter;
alpha ~ normal(0, 10);                        //alpha (intercept)
beta ~ normal(mu, sigma);                     //prior for beta (slope, linked to force of infection)
mu ~ normal(0, 10);                           //prior for global mean
sigma ~ normal(0, 5);                         //prior for variance
gamma ~ normal(0,10);

for(i in 1:N){                                //iterate through patient days (observations)
real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta*colonised[i] + gamma*ecoli_colonized_swab[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
events ~ bernoulli(interval_prob);              //regression
}
generated quantities{
//calculate log likelihood of model - allows for model comparison
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
vector[N_intervals] log_lik;                  //stores log-likelihood of the model
int counter;

for(i in 1:N){                                //iterate through patient days (observations)
real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta*colonised[i] + gamma*ecoli_colonized_swab[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
//calculate log-likelihood of model
for(k in 1:N_intervals){
log_lik[k] = bernoulli_lpmf(events[k]|interval_prob[k]);
}
}'

# MODEL 5
MODEL_5 <-'
data{
int<lower=1> N;                                           //number of observations
int<lower=1> N_intervals;                                 //total number of intervals
int<lower=0, upper=1> events[N_intervals];                //binary outcome 
int<lower=0> colonised[N];                                //patients already colonised in the ward
int<lower=1> N_ward;                                      //number of wards (4)
int<lower=1> ward[N];                                     //references the ward
int<lower=1> interval_index[N];                           //index of intervals
int<lower=1> interval_length[N_intervals];                //interval length
int<lower=0, upper=1> previously_colonised[N];            //binary previously acquired pOXA-48 by other species
}
parameters{
real alpha;                                           //alpha (intercept), now indexed by ward
real beta[N_ward];                                    //beta (coefficient)
real gamma;
real mu;                                              //global mean of intercepts
real<lower=0> sigma;                                  //variance of intercepts
}
model{
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
int counter;
alpha ~ normal(0, 10);                        //alpha (intercept)
beta ~ normal(mu, sigma);                     //prior for beta (slope, linked to force of infection)
mu ~ normal(0, 10);                           //prior for global mean
sigma ~ normal(0, 5);                         //prior for variance
gamma ~ normal(0,10);

for(i in 1:N){                                //iterate through patient days (observations)
real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta[ward[i]]*colonised[i] + gamma*previously_colonised[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
events ~ bernoulli(interval_prob);              //regression
}
generated quantities{
//calculate log likelihood of model - allows for model comparison
vector[N_intervals] interval_prob;            //vector to store interval probabilities
vector[N] p;                                  //vector to store daily probabilities
vector[N_intervals] log_lik;                  //stores log-likelihood of the model
int counter;

for(i in 1:N){                                //iterate through patient days (observations)
real day_odds;
//linear function of covariates, converted onto probability scale
day_odds = exp(alpha + beta[ward[i]]*colonised[i] + gamma*previously_colonised[i]);
p[i] = day_odds / (day_odds + 1);
}

//iterate through intervals, calculate probability per interval
counter = 0;                      
for(k in 1:N_intervals){
vector[(interval_length[k])] interval_temp;     //create a temporary vector of length interval k. Stores 1-p
for(j in 1:(interval_length[k])){
interval_temp[j] = 1-p[(counter+j)];            //fill elements of temp vector with 1-p for each day in interval k
}
//product of each (1-p) in the interval (sum of the logs)
interval_prob[k] = 1-exp(sum(log(interval_temp)));
counter = counter+(interval_length[k]);         //increase the counter by the length of interval k
}
//calculate log-likelihood of model
for(k in 1:N_intervals){
log_lik[k] = bernoulli_lpmf(events[k]|interval_prob[k]);
}
}'


