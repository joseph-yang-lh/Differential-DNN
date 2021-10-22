# Differential-DNN

@ Authors
Tingting Zhou tingting_zhou@berkeley.edu
Joseph Yang joseph_yang@berkeley.edu



Artificial neural networks(ANNs) has been widely used to detect patterns in large dataset. As being a universal function approximators, ANN has also been used in claim pricing with the aim to accelerating the corresponding numerical methods. Most of the work focus on pricing and reducing the MSE of the test dataset. Take vanilla European option for instance, the analytical solution is famous Black-Scholes formula, a 5 layer neural network is able to
deliver a pretty close result to BS solution while with unbalanced MPE along moneyness. 

The project mainly conclude two parts:

In our first stage, we tried to construct a option pricing engine with smooth Greeks,we firstly constructed a 5 layer neural network to price European options, compared to previous work, we trained the model with mesh grid generated data to increase the amount of corner data points and tuned the hyper parameters to generate a more balanced performance along moneyness dimension. Once the pricing model is set up, we extended another ANN model as a weighted smoother, option price and other parameters will be fed into the second ANN model and deliver smooth first order and higher order Greeks. In the training process for the second ANN model, we froze all parameters in the first ANN model and tuned the second ANN model to smooth Greeks until a satisfying result was achieved.

-Folder: DNN_for_Pricing_and_Hedging
--main.ipynb

In our second stage, we introduced stochastic volatility in the underlying process, which has a more complicated approximate analytical pricing solution and exact numerical solving methods. In this part, we introduced SABR stochastic volatility model(Hagan et al.,2002)[1] and its classic Lognormal approximate formula,Floc’h method which is regarded as lognormal expansion to the classic one, finite difference numerical method with boundary constraints and also arbitrage free finite difference method without boundary limitation. We implemented the above four methods and compare them. Finally, we convert price into implied volatility using Black-Scholes lognormal formula, and built ANN model to learn the difference of implied volatility between analyti- cal pricing formula and numerical exact result. Any claim can be easily priced with the closed form pricing plus the difference estimated by well trained ANN model. The pricing efficiency will be highly improved taking advantage of the negligible time spent on analytical calculation and the prediction time cost by ANN.

-Folder:DNN_for_Estimating_IV

--DNN_Model_for_Hagan_Correction.ipynb

--DNN_Model_for_Hagan_Floch_with_Filtered_Data.ipynb

--Folder:fdm_paper_and_matlab_code