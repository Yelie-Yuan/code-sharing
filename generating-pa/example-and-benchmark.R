###################################################
### examples
###################################################
library("wdnet")
args(rpanet)

# set an initial network with 
# two edges (1, 2, 0.5), (3, 4, 2.0)
netwk0 <- list(edgelist = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
               edgeweight = c(0.5, 2.0), directed = TRUE)

# probability of edge creation scenarios
ctr1 <- rpa_control_scenario(alpha = 0.2, beta = 0.6, gamma = 0.2,
                             beta.loop = FALSE, source.first = FALSE)

# weight of new edges
my_rgamma <- function(n) rgamma(n, shape = 5, scale = 0.2)
ctr2 <- ctr1 + rpa_control_edgeweight(sampler = my_rgamma)

# preference functions
ctr3 <- ctr2 + rpa_control_preference(
  ftype = "default",
  sparams = c(1, 2, 0, 0, 1), tparams = c(0, 0, 1, 2, 1))

# generate a PA network
set.seed(12)
netwk3 <- rpanet(nstep = 1e3, initial.network = netwk0, control = ctr3)
names(netwk3)
print(netwk3)

# new edges at each step
ctr4 <- ctr3 + rpa_control_newedge(
  sampler = function(n) rpois(n, 2) + 1,
  snode.replace = FALSE, tnode.replace = FALSE)

# node groups and immediate reciprocal edges
ctr5 <- ctr4 + rpa_control_reciprocal(
  group.prob = c(0.4, 0.6),
  recip.prob = matrix(c(0.4, 0.1, 0.2, 0.5), nrow = 2, byrow = TRUE))

# set an initial network,
# specify group of nodes 1~4 and generate a PA network
netwk0 <- list(
  edgelist = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
  edgeweight = c(0.5, 2), directed = TRUE, nodegroup = c(1, 2, 2, 1)
)
netwk5 <- rpanet(1e3, control = ctr5, initial.network = netwk0)

# set customized preference functions through one-line 
# c++ style expressions
ctr6 <- ctr5 + rpa_control_preference(
  ftype = "customized",
  spref = "log(outs + 1) + 1", tpref = "log(ins + 1) + 1"
)

# set customized preference functions through XPtr
my_spref <- RcppXPtrUtils::cppXPtr(code =
  "double foo(double x, double y) {
    if (x < 1) {
      return 1;
    } else if (x <= 100) {
      return pow(x, 2);
    } else {
      return 200 * (x - 50);
    }
  }")
ctr7 <- rpa_control_preference(
  ftype = "customized", spref = my_spref, tpref = "ins + 1"
)


##########################################################
### code for benchmarks: weighted; default initial network
##########################################################
library(microbenchmark)
library(wdnet)
library(igraph)
library(PAFit)

# rpanet controls
control0 <- rpa_control_scenario(alpha = 1/3,
                                 beta = 1/3,
                                 gamma = 1/3,
                                 xi = 0,
                                 rho = 0) +
  rpa_control_edgeweight(
    sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
  )
control0.5 <- control0 +
  rpa_control_preference(sparams = c(1, 0.5, 0, 0, 0.1),
                         tparams = c(0, 0, 1, 0.5, 0.1))
control1 <- control0 +
  rpa_control_preference(sparams = c(1, 1, 0, 0, 0.1),
                         tparams = c(0, 0, 1, 1, 0.1))
control2 <- control0 +
  rpa_control_preference(sparams = c(1, 2, 0, 0, 0.1),
                         tparams = c(0, 0, 1, 2, 0.1))

# number of replicates
# the 100 replicates were run on HPC; please change it to a smaller number
# if you run the benchmarks on a local machine
times <- 100
# run benchmark
ret_weighted_default_seed <- microbenchmark(
  linear_k05_n1e3 = rpanet(1e3, control = control0.5, method = "linear"),
  linear_k05_n1e4 = rpanet(1e4, control = control0.5, method = "linear"),
  linear_k05_n1e5 = rpanet(1e5, control = control0.5, method = "linear"),
  linear_k05_n1e6 = rpanet(1e6, control = control0.5, method = "linear"),
  linear_k1_n1e3  = rpanet(1e3, control = control1,   method = "linear"),
  linear_k1_n1e4  = rpanet(1e4, control = control1,   method = "linear"),
  linear_k1_n1e5  = rpanet(1e5, control = control1,   method = "linear"),
  linear_k1_n1e6  = rpanet(1e6, control = control1,   method = "linear"),
  linear_k2_n1e3  = rpanet(1e3, control = control2,   method = "linear"),
  linear_k2_n1e4  = rpanet(1e4, control = control2,   method = "linear"),
  linear_k2_n1e5  = rpanet(1e5, control = control2,   method = "linear"),
  linear_k2_n1e6  = rpanet(1e6, control = control2,   method = "linear"),
  linear_k2_n1e7  = rpanet(1e7, control = control2,   method = "linear"),
  binary_k05_n1e3 = rpanet(1e3, control = control0.5, method = "binary"),
  binary_k05_n1e4 = rpanet(1e4, control = control0.5, method = "binary"),
  binary_k05_n1e5 = rpanet(1e5, control = control0.5, method = "binary"),
  binary_k05_n1e6 = rpanet(1e6, control = control0.5, method = "binary"),
  binary_k05_n1e7 = rpanet(1e7, control = control0.5, method = "binary"),
  binary_k1_n1e3 =  rpanet(1e3, control = control1,   method = "binary"),
  binary_k1_n1e4 =  rpanet(1e4, control = control1,   method = "binary"),
  binary_k1_n1e5 =  rpanet(1e5, control = control1,   method = "binary"),
  binary_k1_n1e6 =  rpanet(1e6, control = control1,   method = "binary"),
  binary_k1_n1e7 =  rpanet(1e7, control = control1,   method = "binary"),
  binary_k2_n1e3 =  rpanet(1e3, control = control2,   method = "binary"),
  binary_k2_n1e4 =  rpanet(1e4, control = control2,   method = "binary"),
  binary_k2_n1e5 =  rpanet(1e5, control = control2,   method = "binary"),
  binary_k2_n1e6 =  rpanet(1e6, control = control2,   method = "binary"),
  binary_k2_n1e7 =  rpanet(1e7, control = control2,   method = "binary"),
  times = times
)

##############################################################
### code for benchmarks: weighted; weighted ER initial network
##############################################################

# set initial network: number of nodes 10^4,
# number of edges: 10^4 * 10^4 * 0.01
seed_nnode <- 1e4
p <- 0.01

# generate a new initial network for each iteration of
# all benchmark expressions
i <- 100
new_network <- function() {
  if (i >= 28) {
    i <<- 0
    g <<- erdos.renyi.game(
      seed_nnode, seed_nnode^2 * p, 
      type = "gnm", directed = TRUE, loops = TRUE
    )
    # initial network for wdnet
    edgelist <- as_edgelist(g)
    initial.network <<- list(
      "edgelist" = edgelist, directed = TRUE,
      edgeweight = rgamma(seed_nnode^2 * p, shape = 5, scale = 0.2)
    )
    
    stopifnot(nrow(initial.network$edgelist) == seed_nnode^2 * p)
    stopifnot(max(initial.network$edgelist) == seed_nnode)
  }
  i <<- i + 1
  return(NULL)
}

# number of replicates
# the 100 replicates were run on HPC; please change it to a smaller number
# if you run the benchmarks on a local machine
times <- 100

# rpanet controls
control0 <- rpa_control_scenario(alpha = 1/3, 
                                 beta = 1/3, 
                                 gamma = 1/3, 
                                 xi = 0, rho = 0) + 
  rpa_control_edgeweight(
    sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
  )

control0.5 <- control0 + 
  rpa_control_preference(sparams = c(1, 0.5, 0, 0, 0.1), 
                         tparams = c(0, 0, 1, 0.5, 0.1))
control1 <- control0 + 
  rpa_control_preference(sparams = c(1, 1, 0, 0, 0.1), 
                         tparams = c(0, 0, 1, 1, 0.1))
control2 <- control0 + 
  rpa_control_preference(sparams = c(1, 2, 0, 0, 0.1), 
                         tparams = c(0, 0, 1, 2, 0.1))

# run benchmark
ret_weighted_ER_seed <- microbenchmark(
  linear_k05_n1e3 = rpanet(1e3, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k05_n1e4 = rpanet(1e4, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k05_n1e5 = rpanet(1e5, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k05_n1e6 = rpanet(1e6, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k1_n1e3  = rpanet(1e3, control = control1,   method = "linear", initial.network = initial.network),
  linear_k1_n1e4  = rpanet(1e4, control = control1,   method = "linear", initial.network = initial.network),
  linear_k1_n1e5  = rpanet(1e5, control = control1,   method = "linear", initial.network = initial.network),
  linear_k1_n1e6  = rpanet(1e6, control = control1,   method = "linear", initial.network = initial.network),
  linear_k2_n1e3  = rpanet(1e3, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e4  = rpanet(1e4, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e5  = rpanet(1e5, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e6  = rpanet(1e6, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e7  = rpanet(1e7, control = control2,   method = "linear", initial.network = initial.network),
  binary_k05_n1e3 = rpanet(1e3, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e4 = rpanet(1e4, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e5 = rpanet(1e5, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e6 = rpanet(1e6, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e7 = rpanet(1e7, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k1_n1e3 =  rpanet(1e3, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e4 =  rpanet(1e4, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e5 =  rpanet(1e5, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e6 =  rpanet(1e6, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e7 =  rpanet(1e7, control = control1,   method = "binary", initial.network = initial.network),
  binary_k2_n1e3 =  rpanet(1e3, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e4 =  rpanet(1e4, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e5 =  rpanet(1e5, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e6 =  rpanet(1e6, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e7 =  rpanet(1e7, control = control2,   method = "binary", initial.network = initial.network),
  times = times,
  setup = new_network(),
  control = list(order = "inorder")
)


############################################################
### code for benchmarks: unweighted; default initial network
############################################################

# set an initial network for igraph
g0 <- graph_from_edgelist(matrix(c(1, 2), nrow = 1), directed = TRUE)

# rpanet controls
control0.5 <- rpa_control_scenario(
  alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0) +
  rpa_control_preference(tparams = c(0, 0, 1, 0.5, 0.1))
control1 <- rpa_control_scenario(
  alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0
) +
  rpa_control_preference(tparams = c(0, 0, 1, 1, 0.1))
control2 <- rpa_control_scenario(
  alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0) +
  rpa_control_preference(tparams = c(0, 0, 1, 2, 0.1))

# number of replicates
# the 100 replicates were run on HPC; please change it to a smaller number
# if you run the benchmarks on a local machine
# inaddition, PAFit is time consuming, please consider excluding it
# from the following benchmarks
times <- 100

# run benchmarks
ret_unweighted_default_seed <- microbenchmark(
  PAFit_k05_n1e3 = generate_net(N = 1e3 + 3, alpha = 0.5, num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e3)^0.5 + 0.1),
  PAFit_k05_n1e4 = generate_net(N = 1e4 + 3, alpha = 0.5, num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e4)^0.5 + 0.1),
  PAFit_k05_n1e5 = generate_net(N = 1e5 + 3, alpha = 0.5, num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e5)^0.5 + 0.1),
  PAFit_k1_n1e3 =  generate_net(N = 1e3 + 3, alpha = 1,   num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e3) + 0.1),
  PAFit_k1_n1e4 =  generate_net(N = 1e4 + 3, alpha = 1,   num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e4) + 0.1),
  PAFit_k1_n1e5 =  generate_net(N = 1e5 + 3, alpha = 1,   num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e5) + 0.1),
  PAFit_k2_n1e3 =  generate_net(N = 1e3 + 3, alpha = 2,   num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e3)^2 + 0.1),
  PAFit_k2_n1e4 =  generate_net(N = 1e4 + 3, alpha = 2,   num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e4)^2 + 0.1),
  PAFit_k2_n1e5 =  generate_net(N = 1e5 + 3, alpha = 2,   num_seed = 2, 
                                multiple_node = 1, m = 1, s = 0, offset = 0,
                                custom_PA = c(0:1e5)^2 + 0.1),
  linear_k05_n1e3 = rpanet(1e3, control = control0.5, method = "linear"),
  linear_k05_n1e4 = rpanet(1e4, control = control0.5, method = "linear"),
  linear_k05_n1e5 = rpanet(1e5, control = control0.5, method = "linear"),
  linear_k05_n1e6 = rpanet(1e6, control = control0.5, method = "linear"),
  linear_k1_n1e3  = rpanet(1e3, control = control1,   method = "linear"),
  linear_k1_n1e4  = rpanet(1e4, control = control1,   method = "linear"),
  linear_k1_n1e5  = rpanet(1e5, control = control1,   method = "linear"),
  linear_k1_n1e6  = rpanet(1e6, control = control1,   method = "linear"),
  linear_k2_n1e3  = rpanet(1e3, control = control2,   method = "linear"),
  linear_k2_n1e4  = rpanet(1e4, control = control2,   method = "linear"),
  linear_k2_n1e5  = rpanet(1e5, control = control2,   method = "linear"),
  linear_k2_n1e6  = rpanet(1e6, control = control2,   method = "linear"),
  linear_k2_n1e7  = rpanet(1e7, control = control2,   method = "linear"),
  binary_k05_n1e3 = rpanet(1e3, control = control0.5, method = "binary"),
  binary_k05_n1e4 = rpanet(1e4, control = control0.5, method = "binary"),
  binary_k05_n1e5 = rpanet(1e5, control = control0.5, method = "binary"),
  binary_k05_n1e6 = rpanet(1e6, control = control0.5, method = "binary"),
  binary_k05_n1e7 = rpanet(1e7, control = control0.5, method = "binary"),
  binary_k1_n1e3 =  rpanet(1e3, control = control1,   method = "binary"),
  binary_k1_n1e4 =  rpanet(1e4, control = control1,   method = "binary"),
  binary_k1_n1e5 =  rpanet(1e5, control = control1,   method = "binary"),
  binary_k1_n1e6 =  rpanet(1e6, control = control1,   method = "binary"),
  binary_k1_n1e7 =  rpanet(1e7, control = control1,   method = "binary"),
  binary_k2_n1e3 =  rpanet(1e3, control = control2,   method = "binary"),
  binary_k2_n1e4 =  rpanet(1e4, control = control2,   method = "binary"),
  binary_k2_n1e5 =  rpanet(1e5, control = control2,   method = "binary"),
  binary_k2_n1e6 =  rpanet(1e6, control = control2,   method = "binary"),
  binary_k2_n1e7 =  rpanet(1e7, control = control2,   method = "binary"),
  igraph_k05_n1e3 = sample_pa(1e3 + 2, power = 0.5, zero.appeal = 0.1, start.graph = g0),
  igraph_k05_n1e4 = sample_pa(1e4 + 2, power = 0.5, zero.appeal = 0.1, start.graph = g0),
  igraph_k05_n1e5 = sample_pa(1e5 + 2, power = 0.5, zero.appeal = 0.1, start.graph = g0),
  igraph_k05_n1e6 = sample_pa(1e6 + 2, power = 0.5, zero.appeal = 0.1, start.graph = g0),
  igraph_k05_n1e7 = sample_pa(1e7 + 2, power = 0.5, zero.appeal = 0.1, start.graph = g0),
  igraph_k1_n1e3 =  sample_pa(1e3 + 2, power = 1,   zero.appeal = 0.1, start.graph = g0),
  igraph_k1_n1e4 =  sample_pa(1e4 + 2, power = 1,   zero.appeal = 0.1, start.graph = g0),
  igraph_k1_n1e5 =  sample_pa(1e5 + 2, power = 1,   zero.appeal = 0.1, start.graph = g0),
  igraph_k1_n1e6 =  sample_pa(1e6 + 2, power = 1,   zero.appeal = 0.1, start.graph = g0),
  igraph_k1_n1e7 =  sample_pa(1e7 + 2, power = 1,   zero.appeal = 0.1, start.graph = g0),
  igraph_k2_n1e3 =  sample_pa(1e3 + 2, power = 2,   zero.appeal = 0.1, start.graph = g0),
  igraph_k2_n1e4 =  sample_pa(1e4 + 2, power = 2,   zero.appeal = 0.1, start.graph = g0),
  igraph_k2_n1e5 =  sample_pa(1e5 + 2, power = 2,   zero.appeal = 0.1, start.graph = g0),
  igraph_k2_n1e6 =  sample_pa(1e6 + 2, power = 2,   zero.appeal = 0.1, start.graph = g0),
  igraph_k2_n1e7 =  sample_pa(1e7 + 2, power = 2,   zero.appeal = 0.1, start.graph = g0),
  times = times
)

#######################################################
### code for benchmarks: unweighted; ER initial network
#######################################################

# set initial network: number of nodes 10^4, 
# number of edges 10^4 * 10^4 * 0.01
seed_nnode <- 1e4
p <- 0.01

# generate a new initial network for each iteration of 
# all benchmark expressions
i <- 100
new_network <- function() {
  if (i >= 43) {
    i <<- 0
    # initial network for igraph
    g <<- erdos.renyi.game(seed_nnode, seed_nnode^2 * p, 
                           type = "gnm", directed = TRUE, loops = TRUE)
    
    # initial network for wdnet
    edgelist <- as_edgelist(g)
    initial.network <<- list("edgelist" = edgelist)
    
    stopifnot(nrow(initial.network$edgelist) == seed_nnode^2 * p)
    stopifnot(max(initial.network$edgelist) == seed_nnode)
  }
  i <<- i + 1
  return(NULL)
}

# number of replicates
# the 100 replicates were run on HPC; please change it to a smaller number
# if you run the benchmarks on a local machine
times <- 100

# rpanet controls
control0.5 <- rpa_control_scenario(alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0) +
  rpa_control_preference(tparams = c(0, 0, 1, 0.5, 0.1))
control1 <- rpa_control_scenario(alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0) +
  rpa_control_preference(tparams = c(0, 0, 1, 1, 0.1))
control2 <- rpa_control_scenario(alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0) +
  rpa_control_preference(tparams = c(0, 0, 1, 2, 0.1))

# run benchmark
ret_unweighted_ER_seed <- microbenchmark(
  linear_k05_n1e3 = rpanet(1e3, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k05_n1e4 = rpanet(1e4, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k05_n1e5 = rpanet(1e5, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k05_n1e6 = rpanet(1e6, control = control0.5, method = "linear", initial.network = initial.network),
  linear_k1_n1e3  = rpanet(1e3, control = control1,   method = "linear", initial.network = initial.network),
  linear_k1_n1e4  = rpanet(1e4, control = control1,   method = "linear", initial.network = initial.network),
  linear_k1_n1e5  = rpanet(1e5, control = control1,   method = "linear", initial.network = initial.network),
  linear_k1_n1e6  = rpanet(1e6, control = control1,   method = "linear", initial.network = initial.network),
  linear_k2_n1e3  = rpanet(1e3, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e4  = rpanet(1e4, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e5  = rpanet(1e5, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e6  = rpanet(1e6, control = control2,   method = "linear", initial.network = initial.network),
  linear_k2_n1e7  = rpanet(1e7, control = control2,   method = "linear", initial.network = initial.network),
  binary_k05_n1e3 = rpanet(1e3, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e4 = rpanet(1e4, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e5 = rpanet(1e5, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e6 = rpanet(1e6, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k05_n1e7 = rpanet(1e7, control = control0.5, method = "binary", initial.network = initial.network),
  binary_k1_n1e3 =  rpanet(1e3, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e4 =  rpanet(1e4, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e5 =  rpanet(1e5, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e6 =  rpanet(1e6, control = control1,   method = "binary", initial.network = initial.network),
  binary_k1_n1e7 =  rpanet(1e7, control = control1,   method = "binary", initial.network = initial.network),
  binary_k2_n1e3 =  rpanet(1e3, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e4 =  rpanet(1e4, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e5 =  rpanet(1e5, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e6 =  rpanet(1e6, control = control2,   method = "binary", initial.network = initial.network),
  binary_k2_n1e7 =  rpanet(1e7, control = control2,   method = "binary", initial.network = initial.network),
  igraph_k05_n1e3 = sample_pa(1e3 + seed_nnode, power = 0.5, zero.appeal = 0.1, start.graph = g),
  igraph_k05_n1e4 = sample_pa(1e4 + seed_nnode, power = 0.5, zero.appeal = 0.1, start.graph = g),
  igraph_k05_n1e5 = sample_pa(1e5 + seed_nnode, power = 0.5, zero.appeal = 0.1, start.graph = g),
  igraph_k05_n1e6 = sample_pa(1e6 + seed_nnode, power = 0.5, zero.appeal = 0.1, start.graph = g),
  igraph_k05_n1e7 = sample_pa(1e7 + seed_nnode, power = 0.5, zero.appeal = 0.1, start.graph = g),
  igraph_k1_n1e3 =  sample_pa(1e3 + seed_nnode, power = 1,   zero.appeal = 0.1, start.graph = g),
  igraph_k1_n1e4 =  sample_pa(1e4 + seed_nnode, power = 1,   zero.appeal = 0.1, start.graph = g),
  igraph_k1_n1e5 =  sample_pa(1e5 + seed_nnode, power = 1,   zero.appeal = 0.1, start.graph = g),
  igraph_k1_n1e6 =  sample_pa(1e6 + seed_nnode, power = 1,   zero.appeal = 0.1, start.graph = g),
  igraph_k1_n1e7 =  sample_pa(1e7 + seed_nnode, power = 1,   zero.appeal = 0.1, start.graph = g),
  igraph_k2_n1e3 =  sample_pa(1e3 + seed_nnode, power = 2,   zero.appeal = 0.1, start.graph = g),
  igraph_k2_n1e4 =  sample_pa(1e4 + seed_nnode, power = 2,   zero.appeal = 0.1, start.graph = g),
  igraph_k2_n1e5 =  sample_pa(1e5 + seed_nnode, power = 2,   zero.appeal = 0.1, start.graph = g),
  igraph_k2_n1e6 =  sample_pa(1e6 + seed_nnode, power = 2,   zero.appeal = 0.1, start.graph = g),
  igraph_k2_n1e7 =  sample_pa(1e7 + seed_nnode, power = 2,   zero.appeal = 0.1, start.graph = g),
  times = times,
  setup = new_network(),
  control = list(order = "inorder")
)


##############################################################
### plot results
##############################################################
library(ggplot2)

collect_results <- function(ret) {
  ret <- data.frame("Type" = ret$expr,
                    "Time" = ret$time / 1e9,
                    "Method" = NA,
                    "Power" = NA,
                    "nstep" = NA)
  pb <- txtProgressBar(0, nrow(ret), style = 3)
  for (i in seq_len(nrow(ret))) {
    setTxtProgressBar(pb, i)
    temp <- strsplit(as.character(ret$Type[i]), "_")[[1]]
    ret$Method[i] <- temp[1]
    if (ret$Method[i] == "binary") {
      ret$Method[i] <- "wdnet binary"
    }
    if (ret$Method[i] == "linear") {
      ret$Method[i] <- "wdnet linear"
    }
    
    ret$Power[i] <- switch(temp[2],
                           k05 = "k=0.5",
                           k1 = "k=1",
                           k2 = "k=2")
    ret$nstep[i] <- switch(temp[3],
                           n1e3 = 1e3,
                           n1e4 = 1e4,
                           n1e5 = 1e5,
                           n1e6 = 1e6,
                           n1e7 = 1e7)
    rm(temp)
  }
  aggregate(ret, Time ~ Method + Power + nstep, FUN = median)
}

# unweighted networks
# median runtime for each combination of method, k and nstep
ret_unweighted_default_seed <- collect_results(ret_unweighted_default_seed)
ret_unweighted_default_seed$Type <- "Default Initial Network"
ret_unweighted_ER_seed <- collect_results(ret_unweighted_ER_seed)
ret_unweighted_ER_seed$Type <- "ER Initial Network"

ret_unweighted <- rbind(ret_unweighted_ER_seed, ret_unweighted_default_seed)
ret_unweighted$Type2 <- paste0(ret_unweighted$Method, ret_unweighted$Power)

# plot results
ggplot(ret_unweighted, aes(x = nstep, y = Time, color = Method, group = Method)) +
  geom_path() +
  geom_point(aes(shape = Method), size = 2, alpha = 0.7) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 5),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "bottom",
        legend.text = element_text(margin = margin(r = 0.2, unit = "inch")),
        legend.text.align = 0, 
        aspect.ratio = 1, 
        panel.spacing = unit(0.8, "lines")) +
  facet_grid(Type ~ Power, scales = "free_y") + 
  xlab("Number of Steps") + 
  ylab("Processing Time (Seconds)") + 
  scale_color_manual(values = c("#E69F00", "#009E73", "#000000", "#56B4E9")) +
  scale_shape_manual(values = c(16, 17, 15, 3))

ggsave(filename = "unwegihted_alpha1_median.pdf", width = 8.2, height = 6)


# weighted networks
# get median runtime for each combination of method, k and nstep
ret_weighted_default_seed <- collect_results(ret_weighted_default_seed)
ret_weighted_default_seed$Type <- "Default Initial Network"
ret_weighted_ER_seed <- collect_results(ret_weighted_ER_seed)
ret_weighted_ER_seed$Type <- "ER Initial Network"
ret_weighted <- rbind(ret_weighted_ER_seed, ret_weighted_default_seed)
ret_weighted$Type2 <- paste0(ret_weighted$Method, ret_weighted$Power)

# plot results
ggplot(ret_weighted, 
       aes(x = nstep, y = Time, color = Method, group = Type2)) +
  geom_path() +
  geom_point(aes(shape = Method), size = 2, alpha = 0.7) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 5),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "bottom",
        legend.text = element_text(margin = margin(r = 0.2, unit = "inch")),
        legend.text.align = 0, 
        aspect.ratio = 1, 
        panel.spacing = unit(0.8, "lines")) +
  facet_grid(Type ~ Power, scales = "free_y") + 
  xlab("Number of Steps") + 
  ylab("Processing Time (Seconds)") + 
  scale_color_manual(values = c("#000000", "#56B4E9")) +
  scale_shape_manual(values = c(15, 3))


ggsave(filename = "wegihted_alpha13_median.pdf", width = 8.2, height = 6)


##############################################################
### first 20 nodes dominates the sampling process when k = 2
##############################################################
# generate 1000 networks and collect the percentage of total preference 
# possesed by the first 20 nodes

# weighted networks
ret1 <- c()
ret2 <- c()
k <- 2
control <- rpa_control_preference(tparams = c(0, 0, 1, k, 0.1),
                                  sparams = c(1, k, 0, 0, 0.1)) +
  rpa_control_scenario(alpha = 1/3, beta = 1/3, gamma = 1/3) + 
  rpa_control_edgeweight(
    sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
  )
set.seed(12345)
method <- "linear"
for (i in 1:1000) {
  temp <- rpanet(1e5, control = control, method = method)
  ret1 <- c(ret1, sum(temp$node.attr$spref[1:20]) / 
            sum(temp$node.attr$spref))
  ret2 <- c(ret2, sum(temp$node.attr$tpref[1:20]) / 
            sum(temp$node.attr$tpref))
  rm(temp)
}
fivenum(ret1)
fivenum(ret2)


# unweighted networks
ret <- c()
control <- rpa_control_preference(tparams = c(0, 0, 1, k, 0.1)) +
  rpa_control_scenario(alpha = 1, beta = 0, gamma = 0)
set.seed(12345)
method <- "linear"
for (i in 1:1000) {
  temp <- rpanet(1e5, control = control, method = method)
  ret <- c(ret, sum(temp$node.attr$tpref[1:20]) / 
           sum(temp$node.attr$tpref))
  rm(temp)
}
fivenum(ret)
