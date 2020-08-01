HOW TO USE THIS CODE

Please follow the following steps to compute the MRMI estimator (Gochanour et al., 2020+). Hopefully
this code will soon be implemented in an R package to make things even easier.

Note:

(1) Download files

Download R Code.zip to your computer. Unzip the folder and place all files inside your preferred
R working directory.

(2) Prepare your datafile

Before applying our method, coerce your data into a dataframe with the following conventions:

Covariates, which may be continuous or categorical, should be named/renamed x1, x2, ...., x_n, where n is
the total number of covariates desired.

Treatment variable should be numeric with the values 0 and 1 and named "r"

Your outcome variables for treatments 1 and 0 should be named "y1" and "y0" and should be numeric.

Counsult "Example_data.csv" for an example. This data is a random sample of 500 from the NHANES data used in the
paper, and can be used as a working example to try out the script.

(3) Specify your models

Open the Main.R script. Set working directory to source file location.

The first thing you will need to modify in the script Main.R is the R formula notation for your regression
and propensity score models. Follow the conventions in this example: i.e. "correct" and "incorrect"
regression models should be specified by "ry","rx1" etc. while propensity score models should be specified with
"r", "x1", etc. For real data applications, the truth of the models won't be known, thus in the paper we called
these Set 1 models and Set 2 models based on the set of covariates used.

my_correct_reg=ry~rx1+rx2+rx3+rx4+rx5
my_bad_reg=ry~rx1*rx3+rx1*rx2+rx6+rx7+rx8

my_correct_prop=r~x1+x2+x3+x4+x5
my_bad_prop=r~x1*x3+x1*x2+x6+x7+x8

(4) Edit fMRMI1 and fMRMI0 functions

Edit lines 28-55 and 220-248 of the script Main.R to reflect the number of covariates you have used in the
above models. In the example above, we have 8. Follow the instructions that appear in lines 26-53 of Main.R to do this.

Here is an example of how you would modify lines 220-248 if you have only 5 covariates:

######################

  x1<-dat[,"x1"]
  x2<-dat[,"x2"]
  x3<-dat[,"x3"]
  x4<-dat[,"x4"]
  x5<-dat[,"x5"] # Deleted x6 through x8 lines here
 
  rx1<-x1[r==0]
  rx2<-x2[r==0]
  rx3<-x3[r==0]
  rx4<-x4[r==0]
  rx5=x5[r==0]  # rx6 through rx8 lines here

  
  df_reg=data.frame(ry,rx1,rx2,rx3,rx4,rx5) #removed rx6 through rx8 from this list
  df_prop=data.frame(r,x1,x2,x3,x4,x5) # removed x6 through x8 from this list

  df1=data.frame(rx1=dat$x1,rx2=dat$x2,rx3=dat$x3,rx4=dat$x4,rx5=dat$x5) # remove rx6=dat$x6,....
######################


(5) Run the script

If you'd like to use different values of lambda, K, or L, please specify them in your call to the function FMRMI_final.

Next, run the entire Main.R script from top to bottom, being sure to include the code at the top that reads
in functions from "helper_functions.R" and "MRMI_function.R" 

(6) Examine your results

Your results will be returned to the object 'RES_MRMI' when the code finishes running. Here is a description of
what each part of the results mean:

All are length 9 because we have 9 MRMI estimators.

Numeric:
means: Contains the estimate of the average treatment effect for each of the 9 estimators.
SEs: Standard errors of all such average treatment effects
Low: 95% C.I. for all such estimates of the treatment effect, low value
High: 95% C.I. for the treatment effect, high value
length: Length of all such 95% C.I.s

Logical:
sig: Whether the estimate of the treatment effect is significance, as determined by whether its confidence interval contains zero.







