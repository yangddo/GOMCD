# GOMCD
Gaussian-based oversampling approach for imbalanced and overlapped class adapting the minimum covariance determinant

we propose a Gaussian-based oversampling adapting minimum covariance determinant (GOMCD) to deal with the class imbalance and overlap simultaneously.
By employing distribution-based random number generation, GOMCD is able to generate artificial instances in a way that expands the training sample area, thereby mitigating the risk of overfitting.
GOMCD estimates the distribution by a Gaussian mixture model adapting minimum covariance determinant.
GOMCD conducts re-clustering through outlier detection and removal.
Through re-clustering that mitigates the influence of outliers, GOMCD can limit the expansion of the training sample area.
We defined the degree of class overlap to generate additional instances in the overlapping areas in order to improve the classification of the minority class in those areas. 

GOMCD_parameter.R file calculates the sampling weight for class overlap.
First, estimates the distribution of the minority class by GMM with MCD.
Then, calculates the degree of class overlap to calculate the sampling weight for each component of the GMM.

GOMCD_oversampling.R file generates artificial instances of the minority class.
