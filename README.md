Contact Stress Distribution using Fast Fourier Transform:

A) BaseComplexNumber.java:

Description - Performs basic operations in complex number like adding, subtracting etc.

List of functions:

1) printArray( BaseComplexNumber[] baseComplexNumber )
      Input - An array containing complex number objects baseComplexNumber,
      Description - Prints all the elements in the array baseComplexNumber
      
2) printObject( BaseComplexNumber baseComplexNumber )
      Input - object baseComplexNumber,
      Description - Prints baseComplexNumber
 
3) add( BaseComplexNumber baseComplexNumber )
      Input - object baseComplexNumber,
      Output - object ComplexNum,
      Description - adds an object ComplexNum to baseComplexNumber
      
3) multiply( BaseComplexNumber baseComplexNumber )
      Input - object baseComplexNumber,
      Output - object ComplexNum,
      Description - multiplies an object ComplexNum to baseComplexNumber 
      
4) divide( BaseComplexNumber baseComplexNumber )
      Input - object baseComplexNumber,
      Output - object ComplexNum,
      Description - divides an object ComplexNum to baseComplexNumber  
   
5) multiplyConstant( double constant )
      Input - double,
      Output - object ComplexNum,
      Description - multiples a double constant to baseComplexNumber
    
6) subtract( BaseComplexNumber baseComplexNumber )
      Input - object baseComplexNumber,
      Output - object ComplexNum,
      Description - subtracts an object ComplexNum to baseComplexNumber
      
7) divideInt( int n )
      Input - integer,
      Output - object ComplexNum,
      Description - divides an object ComplexNum with an integer
      
B) ComplexRootUnity

- extends BaseComplexNumber
- one member variable - int n_

Description: argument of complex number is unity. Complex number is expressed in euler form.

List of functions:

getRealPart(), getImaginaryPart(), getMod() - to obtain real part, imaginary part and modulus of complex number 

      
      
C) ComplexModArg

- extends ComplexRootUnity
- one member variable - int r_

Description: argument of complex number is r_. Complex number is expressed in euler form.

List of functions:

getRealPart(), getImaginaryPart() - to obtain real part and imaginary part of complex number


D) ComplexNum

- extends ComplexNum
- two member variables - double real_, double imaginary_

Description: complex number expressed in rectangular form

List of functions:

getRealPart(), getImaginaryPart(), getMod() - to obtain real part, imaginary part and modulus of complex number 

E) FatFourierTransform
Description - Computes coeffients, inverse coefficients and applies fast fourier/inverse fast fourier transform.
List of functions:

1) fftCoeff( int N )
      Input - int N
      Output - BaseComplexNumber[][]
      Description - computes fast fourier transform coefficients. Each element in the array is an object ComplexRootArg
      
2) fftCoeffDomain( double[] x , double[] alpha , int N)
      Input - double[] x , double[] alpha , int N
      Output - BaseComplexNumber[][]
      Description - transalates the coefficients to the domain 
      
3) fftCoeffRealDomain( double[] x , double[] alpha , int N)
      Input - double[] x , double[] alpha , int N
      Output - BaseComplexNumber[][]
      Description - computes fourier cosine terms in the domain
      
4) ifftCoeff( int N )
      Input - int N
      Output - BaseComplexNumber[][]
      Description - computes inverse fourier transform coefficients
      
5) ifftCoeffRealDomain( double[] x , double[] alpha , int N)
      Input - double[] x , double[] alpha , int N
      Output - BaseComplexNumber[][]
      Description -computes inverse fourier cosine terms in the domain
   
6) fft( BaseComplexNumber[] ba )
      Input - BaseComplexNumber[]
      Output - BaseComplexNumber[]
      Description - applies fast fourier transform
      
7) ifft( BaseComplexNumber[] by )
      Input - BaseComplexNumber[]
      Output - BaseComplexNumber[]
      Description - applies inverse fast fourier transform
      
8) lowPass( double[] fs, double[] sa ,double mu0 )
      Input - double[] fs, double[] sa ,double mu0
      Description - low pass filter
           
 E) GuassElimination
      Description - A funtion lsolve(BaseComplexNumber[][] A, BaseComplexNumber[] b) which performs matrix multiplication with elements being complex numbers using GuassElimination method
      
 F) ContactDistribution
      List of functions:
      contactDistribution - Computes contact distribution
      pDistribution - distribution of contact load
      paConverge - convergence of contact distribution
      loadEq - load distribution 
      delta - computation of a paramter
      MeanSqrErr - Mean square error
      sampling - computes sampling points interval in domain
      exT - domain length
      samplingPoints - computes sampling points within domain
      gapFunction
      vpiAlpha - surface displacement profile
      valpha - finds displacement at a point
      addCompMat - adds complex matrices
      MultiplyCompMat - multiplies two complex matrices
      MulColumnMatCompMat - multiplies vector to complex matrix
      MulCompMatDoubColMat - multiplies complex matrix to a double matrix
      trapezoidalRule - numerical method to integrate
      
      
      
      
 
    

    
 
   





