#Code for learning factors using sn-spMF and mapping eQTLs to factors using weighted LR.



sn_spMF/run_MF_train_coph.R: The wrapper for parameter tuning in sn_spMF\
    -- readIn.R: Read in data. Please modify to fit\
    -- Update_FL.R: Implement alternating least squares (ALS) with gradient descent\
        -- fit_F.R: Update factors\
        -- fit_L.R: Update loadings\
        -- compute_obj.R: Compute the objective\
    -- compute_cophenet.R: Compute cophenetic coefficient after multiple runs.\

sn_spMF/run_MF.R: The wrapper for running sn-spMF.\ 
    -- readIn.R: Read in data. Please modify to fit\
    -- Update_FL.R: Implement alternating least squares (ALS) with gradient descent\
	-- fit_F.R: Update factors\
	-- fit_L.R: Update loadings\
	-- compute_obj.R: Compute the objective \


mapping/run_linearReg.R - mapping eQTLs to the factors
