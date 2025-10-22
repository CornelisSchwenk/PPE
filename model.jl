using DataFrames
using DecisionTree

#"Function that trains forest model given parameter matrix X and output Y"
function fr_model(X,Y,ntrees=100,maxdepth=3,mss=10)
	forest_model = RandomForestRegressor(n_trees=ntrees, max_depth = maxdepth,
							     min_samples_split = mss)
	DecisionTree.fit!(forest_model, X,Y)
	return forest_model
end

#"Function to get the impurity based feature score"
function get_impurity(kys,Y)
	X = dp[kys[1]]
	for i in 2:length(kys)
		X = hcat(X,dp[kys[i]])
	end
	fr = fr_model(X,Y)
	for i in 1:length(kys)
		println("$(kys[i])	$(round(impurity_importance(fr)[i],digits=4))")
	end
end

#Function to get the average impurity for 100 model runs
function get_impurity_av(kys,Y)
	X = dp[kys[1]]
	for i in 2:length(kys)
		X = hcat(X,dp[kys[i]])
	end
	
        out = zeros(size(X)[2])
        R2 = 0
        NRMSE = 0
        for i in 1:100
                fr = fr_model(X,Y)
                out = out .+ impurity_importance(fr)
        end
        out = out ./ 100
        
        for i in 1:length(kys)
                println("$(kys[i])      $(round(out[i],digits=4))")
        end

        inds = findall(out .> 0.1)
        println(kys[inds])
        X2 = dp[kys[inds[1]]]
        for i in 2:length(inds)
            X2 = hcat(X2,dp[kys[inds[i]]])
        end

        R2 = 0
        NRMSE = 0
        for i in 1:100
            fr2 = fr_model(X2,Y)
            y_pred = DecisionTree.predict(fr2, X2)
            y_true = Y
            ss_total = sum((y_true .- mean(y_true)) .^ 2)
            ss_residual = sum((y_true - y_pred) .^ 2)
            r_squared = (1 - ss_residual / ss_total)
            R2 = R2 + r_squared
            RMSE = rmsd(y_true,y_pred)
            NRMSE = NRMSE + (RMSE / std(y_true))
        end
        R2 = round(R2 ./ 100,digits=3)
        NRMSE = round(NRMSE ./ 100,digits=3)
        println("R2 = $(R2)")
        println("NRMSE = $(NRMSE)")
end

#"Function to get the R^2 and NRMSE performance of random forest model"
function get_r2(kys,Y)
	X = dp[kys[1]]
	for i in 2:length(kys)
		X = hcat(X,dp[kys[i]])
	end
	fr = fr_model(X,Y)
	y_pred = DecisionTree.predict(fr, X)
	y_true = Y
	ss_total = sum((y_true .- mean(y_true)) .^ 2)
	ss_residual = sum((y_true - y_pred) .^ 2)
	r_squared = round(1 - ss_residual / ss_total,digits=3)
	RMSE = rmsd(y_true,y_pred)
	NRMSE = round(RMSE / std(y_true),digits=3)
	println("R^2 = $(r_squared)")
	println("NRMSE = $(NRMSE)")	
end

