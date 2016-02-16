#                               _             _     _                     _                     _       
#                              | |           | |   (_)                   | |                   | |      
#    _ __ ___     __ _   _ __  | | __   ___  | |_   _   _ __     __ _    | |_    ___     ___   | |  ___ 
#   | '_ ` _ \   / _` | | '__| | |/ /  / _ \ | __| | | | '_ \   / _` |   | __|  / _ \   / _ \  | | / __|
#   | | | | | | | (_| | | |    |   <  |  __/ | |_  | | | | | | | (_| |   | |_  | (_) | | (_) | | | \__ \
#   |_| |_| |_|  \__,_| |_|    |_|\_\  \___|  \__| |_| |_| |_|  \__, |    \__|  \___/   \___/  |_| |___/
#                                                                __/ |                                  
#                                                               |___/                                   

##############################
# R Library: marketingtools  #
##############################

# Developed by: 	Hannes Datta 
					Tilburg University, The Netherlands
					email: h.datta@uvt.nl

# This repository is hosted on https://github.com/hannesdatta/marketingtools.git

############			
# Features #
############

- Computation of Gaussian Copula correction terms, using their ECDF (Park & Gupta 2012)
- Estimation of Feasible Generalized Least Squares (FGLS) with unbalanced panels
- Estimation of market share attraction models (Fok 2001)
  o Transformation of wide panels into their base-brand representation
  o Support for homogenous and heterogenous coefficients
  o Restricted competition; fully extended model not implemented
  o Estimation via itersur() (FGLS)
  
