ccol=~/src/prog/cCOLUMBUS
echo molint
${ccol}/molint < int.in > int.out
echo cdip
${ccol}/cdip < int.in > cdip.out
<<<<<<< HEAD
echo pkfs
${ccol}/pkfl < dum.in > pkfl.out
echo cscf
${ccol}/cscf < scf.in > scf.out
echo CISD_ground
${ccol}/molci < ci1.in > ci1.out
mv CSF CSFG
mv CIVEC CIVECG
echo CISD_EE
=======
echo pkfl
${ccol}/pkfl < dum.in > pkfl.out
echo scf
${ccol}/cscf < scf.in > scf.out

echo "CI(Groud)"
${ccol}/molci < ci1.in > ci1.out
mv CSF CSFG
mv CIVEC CIVECG
echo "CI(EE)"
>>>>>>> ac9f05299e006c468a7707e24245b07a2de4b2d2
${ccol}/molci < ci2.in > ci2.out
