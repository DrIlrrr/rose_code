pwd
echo ' '
echo 'on the machine:'
hostname


#matlab -nodisplay -r "run_cain2, quit" > report.out
#/usr/local/MATLAB/R2013a/bin/matlab -nodisplay -r "run_cain2, quit" > report.out
#/usr/local/MATLAB/R2013a/bin/matlab  -nojvm -r run_cain2 -logfile report.out </dev/null
/usr/local/MATLAB/R2016b/bin/matlab  -nojvm -r run_using_elegant_beam_cain2 -logfile report.out </dev/null 
