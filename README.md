# I4I5_model

### Instruction to use I4I5 model

The umat was implemented based on the I4I5 model [1-4]

1. In the .inp file material section, specify the input as follows:

   *Material, name=FengUmat

   *User Material, Constants=7

   13.0, 50, 1.0, 1.0e3, 0, 1, 0.0,


2. From left to right, the 7 constants are:

   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\Psi&space;=&space;\frac{\mu}{2}[(I_1-3)&plus;\zeta(I_4-1)^2&plus;\phi&space;I_5^*]" title="\Psi = \frac{\mu}{2}[(I_1-3)+\zeta(I_4-1)^2+\phi I_5^*]" />
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;\mu" title="\inline \mu" />- shear modulus
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;\zeta" title="\inline \zeta" />– stretch parameter
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;\phi" title="\inline \phi" /> – shear parameter
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;D" title="\inline D" /> – for implementing the incompressibility
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;x_A" title="\inline x_A" />– x component of the fiber vector
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;y_A" title="\inline y_A" />– y component of the fiber vector
   
   <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\inline&space;z_A" title="\inline z_A" /> – z component of the fiber vector



3. Specify the Gaussian point in Umat

   Check line 12 and 13, RTOTAL and DROTINC, the first number of the matrix has to be set to the total number of Gaussian points in your simulation. In this example, 8 points were set.



### References

1. Feng, Y., et al., *Characterizing white matter tissue in large strain via asymmetric indentation and inverse finite element modeling.* Journal of the Mechanical Behavior of Biomedical Materials, 2017. **65**: p. 490-501.

2. Feng, Y., et al., *On the accuracy and fitting of transversely isotropic material models.* Journal of the Mechanical Behavior of Biomedical Materials, 2016. **61**: p. 554-566.

3. Feng, Y., et al., *Measurements of mechanical anisotropy in brain tissue and implications for transversely isotropic material models of white matter.* Journal of the Mechanical Behavior of Biomedical Materials, 2013. **23**: p. 117-132.

4. Feng, Y., et al., *A computational study of invariant I5 in a nearly incompressible transversely isotropic model for white matter.* Journal of Biomechanics, 2017. **57**: p. 146-151.

 