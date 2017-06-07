pwd
echo ' '
echo 'on the machine:'
hostname

cd /home/camilla/Desktop/other/new_code/code_for_cardarelli_new_highten_5.8_MeV
pwd
/usr/local/MATLAB/R2013a/bin/matlab -nodisplay -nojvm -r "photons_plots_not_function, all_in_one, quit"

mail dr.ilrrr@gmail.com <<-EOT
this job for 5.8 MeV
is finished
EOT