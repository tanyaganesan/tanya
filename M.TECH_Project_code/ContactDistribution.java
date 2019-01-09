import java.util.ArrayList;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;


public class ContactDistribution 
{
  public static BaseComplexNumber[] contactDistribution( double R, double b, double P, double a, double delta, double E1, double mu1, double E2, double mu2, int N )
  {
    double[][] samplePoints = samplingPoints( N , a ) ;
    double exTa = exT( a ) ;
    System.out.println( "exTa" + exTa ) ;
    double bandwidth = sampling( a , N ) ;
    System.out.println( "bandwidth" + bandwidth ) ;
    double dx = 1 / ( 2 * bandwidth ) ;
    System.out.println( "dx" + dx ) ;
    double dAlpha = 1 / ( 2 * exTa ) ;
    System.out.println( "dAlpha" + dAlpha ) ;
    ////Print out sample points 
//    for( int i = 0 ; i < 2 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        System.out.println( "samplePoints[i][j]]" + samplePoints[i][j] ) ;
//      }
//    }
    /////Surface displacement in layer of thickness 2b. Origin at middle of the layer.
    double[] surfaceDisp1 = new double[N+1] ;
    surfaceDisp1 = vpiAlpha( mu1, E1 , b , N , a ) ;
    
    //Print out surface displacement in the layer.
//    for( int i = 0 ; i < N + 1 ; i++ )
//    {
//      System.out.println( "surfaceDisp1[i]" + surfaceDisp1[i] ) ;
//    }
//    
    /////Surface displacement of the cylindrical indenter
    double[] surfaceDisp2 = new double[N+1] ;
    double[] surfaceDisp2temp = new double[N+1] ;
    surfaceDisp2temp = vpiAlpha( mu2, E2 , b , N , a ) ;
    for( int i = 0 ; i < N+1 ; i++ )
    {
      surfaceDisp2[i] = - surfaceDisp2temp[i] ;
    } 
    
    //Print out surface displacement in the indenter.
//    for( int i = 0 ; i < N + 1 ; i++ )
//    {
//      System.out.println( "surfaceDisp2[i]" + surfaceDisp2[i] ) ;
//    }
    
    
    //////Displacement at y = -b 
   
//    double[] vy = new double[N+1] ;
//    vy = valpha( mu2, E2 , b , N , a ) ;
    
    ///////For any type of p(x) loading 
    //BaseComplexNumber[][] fftcoeffN = FastFourierTransform.fftCoeffDomain( samplePoints[0] , samplePoints[1] , N ) ;
    
//    /////// Multiply each element with dx 
//    for( int i = 0 ; i < N + 1 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        fftcoeffN[i][j] = fftcoeffN[i][j].multiplyConstant( dx ) ;
//      }    
//    }
    
    //BaseComplexNumber[][] ifftcoeffN = FastFourierTransform.ifftCoeffDomain( samplePoints[0] , samplePoints[1] , N ) ;
    
//    /////// Multiply each element with dAlpha 
//    for( int i = 0 ; i < N + 1 ; i++ )
//    {
//      for( int j = 0 ; j < N + 1 ; j++ )
//      {
//        ifftcoeffN[i][j] = ifftcoeffN[i][j].multiplyConstant( dAlpha ) ;
//      }      
//    }
    
    //////Considering p(x) to be real and even, only Fourier cosine coefficients are considered 
    BaseComplexNumber[][] fftcoeffN = FastFourierTransform.fftCoeffRealDomain( samplePoints[0] , samplePoints[1] , N ) ;
    
    /////// Multiply each element with dx 
    for( int i = 0 ; i < N + 1 ; i++ )
    {
      for( int j = 0 ; j < N+1 ; j++ )
      {
        fftcoeffN[i][j] = fftcoeffN[i][j].multiplyConstant( dx ) ;
        //System.out.println( " fftcoeffN " +  fftcoeffN[i][j] ) ;
      }    
    }
    
    BaseComplexNumber[][] ifftcoeffN = FastFourierTransform.ifftCoeffRealDomain( samplePoints[0] , samplePoints[1] , N ) ;
    
    /////// Multiply each element with dAlpha 
    for( int i = 0 ; i < N + 1 ; i++ )
    {
      for( int j = 0 ; j < N + 1 ; j++ )
      {
        ifftcoeffN[i][j] = ifftcoeffN[i][j].multiplyConstant( dAlpha ) ;
        //System.out.println( " ifftcoeffN " +  ifftcoeffN[i][j] ) ;
      }      
    }
    
    //////[T]= [ [F]^-1 ]{ vpi(alpha) }[F] . vpi column matrix multiplied to each column of Fourier coefficient matrix.
    BaseComplexNumber[][] T1temp = MulColumnMatCompMat( surfaceDisp1 , fftcoeffN ) ;
    BaseComplexNumber[][] T1 = MultiplyCompMat( ifftcoeffN , T1temp ) ;
    
    //Print out T1 matrix
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        System.out.println( "T1[i][j]]" + T1[i][j] ) ;
//      }
//    }
    
    
    BaseComplexNumber[][] T2temp = MulColumnMatCompMat( surfaceDisp2 , fftcoeffN ) ;
    BaseComplexNumber[][] T2 = MultiplyCompMat( ifftcoeffN , T2temp ) ;
    
    //Print out T2 matrix
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        System.out.println( "T2[i][j]]" + T2[i][j] ) ;
//      }
//    }
    
    BaseComplexNumber[][] T = addCompMat( T1 , T2 ) ;
    //Print out T matrix
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        System.out.println( "T[i][j]]" + T[i][j] ) ;
//      }
//    }
    
    ////Gap function 
    BaseComplexNumber[] ha = new BaseComplexNumber[N+1] ;
    ha = gapFunction( R , delta , a , N ) ;
    
    //Print out gap function
    for( int i = 0 ; i < N + 1 ; i++ )
    {
      System.out.println( "ha[i]" + ha[i] ) ;
    }
    
    /////// { p(x) } = [ [T]^-1 ]{h(x) - delta } Solved using Gauss elimination method 
    
    /////// p(x) is even function. Only half of the equations are required.
    BaseComplexNumber[][] Thalf = new BaseComplexNumber[ ( N+1 ) / 2 ][ ( N+1 ) / 2 ] ;
    for( int i = 0 ; i < ( N+1 ) / 2  ; i++ )
    {
      for( int j = 0 ; j < ( N+1 ) / 2 ; j++ )
      {
        Thalf[i][j] = T[i][j] ;
      }
    }
    
    BaseComplexNumber[] hHalf = new BaseComplexNumber[ ( N+1 ) / 2 ] ;
    for( int i = 0 ; i < ( N+1 ) / 2 ; i++ )
    {
      hHalf[i] = ha[i] ;
    }
    BaseComplexNumber[] pa = new BaseComplexNumber[ ( N+1 ) / 2] ;
    pa = GaussElimination.lsolve( Thalf , hHalf ) ;
    
//    BaseComplexNumber[] pa_tensile = new BaseComplexNumber[ ( N+1 ) / 2] ;
//    for( int i = 0 ; i < (N + 1) / 2 ; i++ )
//    {
//      if( pa[i].getRealPart() < 0.0 )
//      {
//        pa_tensile[i] = new ComplexNum( 0.0 , 0.0 ) ;
//      }
//      else
//      {
//        pa_tensile[i] = pa[i] ;
//      }
//    }
    return pa ;
    //return pa_tensile ;
  }
  
  ///////
  public static void pDistribution( double R , double b , double P , double a , double delta , double E1 , double mu1 , double E2 , double mu2 , int  N )
  {
    ///// convergence of contact length 
    double aModified = a ;
    //BaseComplexNumber[] pa = new BaseComplexNumber[ (N + 1) / 2 ] ;
    for( int i = 0 ; i < 15 ; i++ )
    {
     //System.out.println( " seen " ) ;
     double bandwidth = ( N / ( 2.0 * aModified ) ) ;
     BaseComplexNumber[] pa = contactDistribution( R , b , P , aModified , delta , E1 , mu1 , E2 , mu2 , N ) ;
     for( int j = 0 ; j < pa.length ; j++ )
     {
       System.out.println( "pa[j]" + pa[j] ) ;
     }
     aModified = loadEq( pa , bandwidth , P , aModified ) ;
     System.out.println( "aModified" + aModified ) ;
    }
  
    
//    for( int i = 0 ; i < pa.length ; i++ )
//    {
//      System.out.println( "pa[i]" + pa[i] ) ;
//    }
  }
  
  public static void paConverge( double R, double b, double P, double a, double delta, double E1, double mu1, double E2, double mu2 )
  {
    ///// pa is purely real
    ArrayList<double[]> pList = new ArrayList< double [] >() ;
    for( int n = 9 ; n < 1000 ; n = 2*n )
    {
      BaseComplexNumber[] pa = contactDistribution( R , b , P , a , delta, E1 , mu1 , E2 , mu2 , n ) ;
      double[] p = new double[ pa.length ]  ;
      for( int i = 0 ; i < pa.length ; i++ )
      {
        p[i] = pa[i].getRealPart() ;
      }
      pList.add( p ) ;
    }
    
    double[] mse = new double[pList.size()];
    for( int j = 0 ; j < pList.size() -  1 ; j++ )
    {
     mse[j]=  MeanSqrErr(  pList.get( j ) ,   pList.get( j + 1 ) )  ;
     System.out.println( " mse[j] " + mse[j]) ;
    }
    
  } 
  
  public static double loadEq( BaseComplexNumber[] pa , double bandwidth , double P , double a )
  {
    ////// P is proportional to a^2.
    double[] p = new double[ pa.length ]  ;
    for( int i = 0 ; i < pa.length ; i++ )
    {
      p[i] = pa[i].getRealPart() ;
    }
    double I = trapezoidalRule( p , bandwidth ) ;
    System.out.println( " I "+ I ) ;
    double aModified = Math.sqrt( P / ( I / Math.pow(  a , 2.0 ))  ) ;
    return aModified ;
  }
  
  public static double[] delta( double R, double delta, double a, int N , double[] surfaceDisp1 , double[] surfaceDisp2 )
  {
    BaseComplexNumber[] ha = gapFunction( R , delta , a , N ) ;
    double[] h = new double[ha.length] ;
    for( int i = 0 ; i < ha.length ; i++ )
    {
     h[i] = ha[i].getRealPart() ; 
    }
    double[] hNew = new double[ha.length] ;
    for( int i = 0 ; i < ha.length ; i++ )
    {
      hNew[i] = hNew[i] + delta ; 
    }
    double[] deltaNew = new double[ha.length] ;
    for( int i = 0 ; i < ha.length ; i++ )
    {
      deltaNew[i] = hNew[i] - surfaceDisp1[i] - surfaceDisp2[i] ; 
    }
    return deltaNew ;
  }
  
  //////Computes MSE between two doubles of length N and 2N.
  public static double MeanSqrErr( double[] a1 , double[] a2 )
  {
    /////////// a1 is of length N. a2 is of length 2N.
    double error = 0.0 ;
    int j = 0 ;
    for( int i = 0 ; i < a1.length ; i++ )
    {
      double temp = a1[i] - a2[j] ;
      error = ( error + Math.pow( temp , 2.0 ) ) / Math.pow( 10 , 90 )  ;
      j = j + 2 ;
    }
    double mse = ( error / a1.length ) ;
    System.out.println( "error" + error ) ;
    return mse ;
  }
  
  public static double sampling( double a, int N )
  {
     ////// N = 2a/ ( 1 / 2B ) 
     double exTa = exT( a );
     double bandwidth = ( N / ( 2.0 * exTa ) ) ;
     return bandwidth ;
  }
  
  public static double exT( double a )
  {
    double exT = 1.0 ;
    double exTa = exT * a ;
    return exTa ;
  }
  
  //////Finds sample points in both frequency and x domain
  ////// Given N should be odd
  ///////N + 1 sampling points will be developed which is even. Half of the points to the left of x = 0 and half to the right.
  public static double[][] samplingPoints( int N, double a )
  {
    double bandwidth = sampling(  a, N ) ;
    double exTa = exT( a ) ;
    //System.out.println( bandwidth );
    double[][] samplingPoints = new double[2][N+1] ;
    //frequency domain from -B to B
    samplingPoints[0][0] = -(bandwidth/2) ;
    for( int i = 1 ; i < N + 1  ; i++ )
    {
      samplingPoints[0][i] = -(bandwidth/2) + ( i ) / ( 2*exTa );
    }
    
    //x domain from -a to a 
    samplingPoints[1][0] = -exTa ;
    for( int i = 1 ; i < N + 1 ; i++ )
    {
      samplingPoints[1][i] = -exTa +( i ) / ( bandwidth );
      System.out.println( "samplingPoints[1][i]" + samplingPoints[1][i] ) ;
    }
    return samplingPoints ;
  }
  
  ////// Gap function is assumed to be a second order function. 0 outside contact length.
  public static BaseComplexNumber[] gapFunction( double R, double delta, double a, int N )
  {
    double[][] samplingPoints = new double[2][N+1] ;

    double[] samplingPointsX = new double[N+1] ;
    samplingPoints = samplingPoints( N , a ) ;
    for( int i = 0 ; i < N+1  ; i++ )
    {
      samplingPointsX[i] = samplingPoints[1][i] ;
      //System.out.println( " samplingPointsX[i] " + samplingPointsX[i]  ) ;
    }    
    double[] ha = new double[ N+1 ] ;
    for( int j = 0 ; j < N + 1 ; j++ )
    {
      if( samplingPointsX[j] < a && samplingPointsX[j] > -a )
      {
        ha[j] = ( 1.0 / R ) * Math.pow( samplingPointsX[j] , 2.0 ) -  delta ;
      }
      else
      {
        ha[j] = 0.0 ;
      }
    }
    BaseComplexNumber[] haComplex = new BaseComplexNumber[N+1] ;
    for( int i = 0 ; i < N+1 ; i++)
    {
      ComplexNum temp = new ComplexNum( ha[i] , 0.0 ) ; 
      haComplex[i] = temp ;
    }
    return haComplex ;
  }
  
  ///////// Surface displacements as functions of alpha ( inside integral) 
  
  //// vpialpha is an even function
  public static double[] vpiAlpha( double mu, double E, double b, int N, double a )
  {
    double[][] samplingPoints = new double[2][ N + 1 ] ;
    double[] sampPtsA = new double[N + 1] ;
    samplingPoints = samplingPoints( N , a ) ;
    for( int i = 0 ; i < N + 1; i++ )
    {
      ////////Frequency domain
      sampPtsA[i] = samplingPoints[0][i] ;
    }
    
    double[] vpiAlphaA = new double[ N + 1 ] ;
    for( int j = 0 ; j < N + 1 ; j++ )
    {
      ///// when alpha = 0 . Limit = b/2 
      if( sampPtsA[j] == 0 )
      {
        double constant = ( ( 1.0 + mu )/ ( 2.0 * Math.PI * E )) * ( 2.0 * ( 1.0 - mu ) ) ;
        vpiAlphaA[j] = constant* b / 2 ;
      }
      else
      {
        double temp = sampPtsA[j] * b ;
        double sinhTemp = Math.sinh( temp ) ;
        double coshTemp = Math.cosh( temp ) ;
        double constant = ( ( 1.0 + mu )/ ( 2.0 * Math.PI * E )) * ( 2.0 * ( 1.0 - mu ) ) ;
        Double nr = (- 1.0 /sampPtsA[j] ) * Math.pow( sinhTemp , 2 ) ;
        //System.out.println( " nr " +  nr ) ;
        Double dr = temp + coshTemp * sinhTemp ;
        //System.out.println( " dr " +  dr ) ;
        boolean isNrNan = nr.isNaN() ;
        boolean isDrNan = dr.isNaN() ;
        if( isNrNan &&  isDrNan || nr == Double.POSITIVE_INFINITY && dr == Double.NEGATIVE_INFINITY || nr == Double.NEGATIVE_INFINITY  && dr == Double.POSITIVE_INFINITY )
        {
          vpiAlphaA[j] = constant ;
        }
        else
        {
          vpiAlphaA[j] = constant * nr / dr  ; 
        }       
      }
    }
    return vpiAlphaA ;
  }
  
  ///////////Finds displacement at y = -b 
  public static double[] valpha( double mu, double E, double b, int N, double a )
  {
    double[][] samplingPoints = new double[2][ N + 1 ] ;
    double[] sampPtsA = new double[N + 1] ;
    samplingPoints = samplingPoints( N , a ) ;
    for( int i = 0 ; i < N + 1; i++ )
    {
      sampPtsA[i] = samplingPoints[0][i] ;
      //System.out.println("sampPtsA[i]" + sampPtsA[i] ) ;
    }
    
    double[] valpha = new double[ N + 1 ] ;
    for( int j = 0 ; j < N + 1 ; j++ )
    {
      ///// when alpha = 0 . Limit = b/2 
      if( sampPtsA[j] == 0 )
      {
//        double constant = ( ( 1.0 + mu )/ ( 2.0 * Math.PI * E )) * ( 2.0 * ( 1.0 - mu ) ) ;
//        vpiAlphaA[j] = constant* b / 2 ;
        valpha [j] = 0.0 ;
      }
      else
      {
          double BNum = ( sampPtsA[j] * Math.sinh( sampPtsA[j]*b ) ) ;
          double BDen =  ( 2.0 * Math.pow( sampPtsA[j] , 2.0 ) )*( sampPtsA[j]*b + Math.cosh( sampPtsA[j]*b )*Math.sinh( sampPtsA[j]*b )) ;
          double B = BNum / BDen ;
          double A = (  ( 1 /  sampPtsA[j] ) + ( b * Math.cosh( sampPtsA[j]*b )/ Math.sinh( sampPtsA[j]*b )) )*B  ;
          double c1 = ( ( 1.0 + mu )/ ( 2.0 * Math.PI * E )) ;
          double c2 = 1 - mu ;
          double c3 = 2 + mu ;
          double term11 = ( - Math.pow( sampPtsA[j] , 3.0 ) * ( A - b*B ) * Math.pow( Math.E , sampPtsA[j]*b ) ) ;
          double term12 = Math.pow( sampPtsA[j] , 2.0  ) * B * Math.pow( Math.E , sampPtsA[j]*b ) ;
          double term13 = -sampPtsA[j]*( 1 - sampPtsA[j] ) * B * Math.pow( Math.E , sampPtsA[j]*b ) ;
          double term1 = c2 * ( term11 + term12 + term13 ) ;
          double term21 = ( A - b* B )* Math.pow( Math.E , sampPtsA[j]*b ) ;
          double term22 = ( 1 - sampPtsA[j] ) * B * Math.pow( Math.E , sampPtsA[j]*b ) ;
          double term2 = c3 * Math.pow( sampPtsA[j] , 2.0 ) * ( term21 + term22  ) ;
          valpha [j] = c1 * ( term1 + term2 ) / Math.pow( sampPtsA[j] , 2.0 )  ;
      }
    }
    return valpha ;
  }
  
  //////// adds two matrices. Matrix elements being complex.
  public static BaseComplexNumber[][] addCompMat( BaseComplexNumber[][] mat1 , BaseComplexNumber[][] mat2 )
  {  
    int row = mat1.length ;
    int column = mat1[0].length ;
    BaseComplexNumber[][] addedMat = new BaseComplexNumber[row][column] ;
    for( int i = 0 ; i < row ; i++ )
    {
      for( int j = 0 ; j < column ; j++ )
      {
        addedMat[i][j] = mat1[i][j] ;
        addedMat[i][j] = addedMat[i][j].add( mat2[i][j] ) ;
      }
    }
    return addedMat ;
  }
  
  //////// multiplies two matrices. Matrix elements being complex.
  public static BaseComplexNumber[][] MultiplyCompMat( BaseComplexNumber[][] mat1 , BaseComplexNumber[][] mat2 )
  {  
    int row = mat1.length ;
    int column = mat2[0].length ;
    BaseComplexNumber[][] MultiplyMat = new BaseComplexNumber[row][column] ;
    for( int i = 0 ; i < row ; i++ )
    {
      for( int j = 0 ; j < column ; j++ )
      {
        MultiplyMat[i][j] = mat1[i][j] ;
        MultiplyMat[i][j] = MultiplyMat[i][j].multiply( mat2[i][j] ) ;
      }
    }
    return MultiplyMat ;
  }
  
////////multiplies a double column matrix to a each column of a complex matrix.
  public static BaseComplexNumber[][] MulColumnMatCompMat( double[] mat1 , BaseComplexNumber[][] mat2 )
  {  
    int row = mat1.length ;
    int column = mat2[0].length ;
    BaseComplexNumber[][] MultiplyMat = new BaseComplexNumber[row][column] ;
    for( int j = 0 ; j < column ; j++ )
    {   
      for( int i = 0 ; i < row ; i++ )
      {
        MultiplyMat[i][j] = mat2[i][j] ;
        MultiplyMat[i][j] = MultiplyMat[i][j].multiplyConstant( mat1[i] ) ;
      }
    }
    return MultiplyMat ;
  }
  
  //////////multiplies a complex matrix to a double column matrix.
  public static BaseComplexNumber[][] MulCompMatDoubColMat( BaseComplexNumber[][] mat1 , double[] mat2 )
  {  
    int row = mat1.length ;
    int column = mat2.length ;
    BaseComplexNumber[][] MultiplyMatColMat = new BaseComplexNumber[row][column] ;
    BaseComplexNumber temp = new ComplexNum( 0.0 , 0.0 ) ;
    for( int j = 0 ; j < column ; j++ )
    {   
      for( int i = 0 ; i < row ; i++ )
      {
        for( int k = 0 ; k < column ; k++ )
        {
           temp = mat1[i][j].multiplyConstant( mat2[j] ) ;
           MultiplyMatColMat[i][k] = MultiplyMatColMat[i][k].add( temp ) ;
        }
      }
    }
    return MultiplyMatColMat ;
  }
  
  ///////// Finds integral using trapezoidal rule given a function at discrete points.
  public static double trapezoidalRule( double[] pa , double bandwidth )
  {
    double h =  1.0/bandwidth ;
    //System.out.println("h" + h ) ;
    double n = pa.length ;
    //System.out.println("n" + n ) ;
    double Ix0x2 = 0.0 ;
    for( int i = 0 ; i < n - 2 ; i++ )
    {
     double Ixi = h* ( ( pa[i+1] + pa[i] ) / 2.0 ) ;
     //System.out.println("Ixi" + Ixi ) ;
     Ix0x2 = Ix0x2 + Ixi ;
     //System.out.println("Ix0x2" + Ix0x2 ) ;
    }
    //even function hence multiplied by 2
    double I = 2.0*Ix0x2 ;
    return I  ;   
  }
  
  public static void main( String[] arg )
  {
    // Testing
    int N = 1023 ;
    double a = 1.3*Math.pow(  10.0 , -8.0 ) ;
    //double a = Math.pow(  10.0 , -4.0 ) ;
    double delta = 9.0*Math.pow(  10.0 , -9.0 ) ;
    double R = 1.5*Math.pow(  10.0 , -3.0 ) ;
    double b = 1.0*Math.pow(  10.0 , -4.0 ) ; ////// 2b is the thickness of the layer 
    double E1 = 2.42*Math.pow(  10.0 , 11.0 ) ; ////// Young's modulus of layer 
    double mu1 = 0.2963 ;
    double P = 5.0 ;
    double E2 = 210.0*Math.pow(  10.0 , 9.0 ) ; ////// Young's modulus of of indenter 
    double mu2 = 0.3 ;
//    
    BaseComplexNumber[] pa = new BaseComplexNumber[N+1] ;
    
    pa = contactDistribution( R , b , P , a , delta , E1 , mu1 , E2 , mu2 , N ) ;
    for( int i = 0 ; i < ( N+1 )/2  ; i++ )
    {
      System.out.println("i" + i + "pa" + pa[i]) ;
    }
    
    
    
//    BaseComplexNumber[] haTest = new BaseComplexNumber[ N + 1 ] ;
//    haTest = gapFunction( R, delta, a , N ) ;
//    for( int i = 0 ; i < haTest.length ; i++)
//    {
//      System.out.println( haTest[i] ) ;
//    }
////  int N = 8 ;
//
////  double a = Math.pow(  10.0 , -4.0 ) ;
//
//
//  double bandwidth = sampling( a , N ) ;
//  System.out.println( "bandwidth" + bandwidth ) ;
//  double[][] testSamplePts = new double[2][N+1] ;
//  testSamplePts = samplingPoints( N , a ) ;
//  for( int i = 0 ; i < 2 ; i++ )
//  {
//    for( int j = 0 ; j < N + 1  ; j++ )
//    {
//      System.out.println( "i" + i );
//      System.out.println( "j" + j );
//      System.out.println( "testSamplePts[i][j] " +  testSamplePts[i][j] );
//    }
//  }
////  int N = 8 ;
//  double[] viAlphaTest = new double[ N + 1] ;
//  viAlphaTest = vpiAlpha( mu1, E1, b, N, a ) ;
//  for( int i = 0 ; i < viAlphaTest.length ; i++)
//  {
//    System.out.println( viAlphaTest[i] ) ;
//  }
//  
////  BaseComplexNumber[][] ifftcoeffN = FastFourierTransform.ifftCoeffDomain( testSamplePts[0] , testSamplePts[1] , N ) ;
////for( int i = 0 ; i < ifftcoeffN.length ; i++)
////{
////  for( int j = 0 ; j < ifftcoeffN[0].length ; j++ )
////  {
////    System.out.println("i" + i ) ;
////    System.out.println("j" + j ) ;
////    System.out.println( ifftcoeffN[i][j] ) ;
////  }
////  
////}
//  
////    BaseComplexNumber[][] mat1Test = new BaseComplexNumber[2][2] ;
////    ComplexNum temp =  new ComplexNum( 1.0, 1.0 ) ;
////    mat1Test[0][0] = temp ;
////    mat1Test[0][1] = temp ;
////    mat1Test[1][0] = temp ;
////    mat1Test[1][1] = temp ;
////    BaseComplexNumber[][] mat2Test = new BaseComplexNumber[2][2] ;
////    mat2Test[0][0] = temp ;
////    mat2Test[0][1] = temp ;
////    mat2Test[1][0] = temp ;
////    mat2Test[1][1] = temp ;
////    BaseComplexNumber[][] addCompMatTest = new BaseComplexNumber[2][2] ;
////    addCompMatTest = addCompMat( mat1Test, mat2Test ) ;
////    for( int i = 0 ; i < 2 ; i ++ )
////    {
////      for( int j = 0 ; j < 2 ; j++ )
////      {
////        BaseComplexNumber.printObject( addCompMatTest[i][j] ) ;
////      }      
////    }
////    BaseComplexNumber[][] mat1Test = new BaseComplexNumber[2][2] ;
////    ComplexNum temp =  new ComplexNum( 1.0, 0.0 ) ;
////    mat1Test[0][0] = temp ;
////    mat1Test[0][1] = temp ;
////    mat1Test[1][0] = temp ;
////    mat1Test[1][1] = temp ;
////    BaseComplexNumber[][] mat2Test = new BaseComplexNumber[2][2] ;
////    mat2Test[0][0] = temp ;
////    mat2Test[0][1] = temp ;
////    mat2Test[1][0] = temp ;
////    mat2Test[1][1] = temp ;
////    BaseComplexNumber[][] MultiplyCompMatTest = new BaseComplexNumber[2][2] ;
////    MultiplyCompMatTest = MultiplyCompMat( mat1Test, mat2Test ) ;
////    for( int i = 0 ; i < 2 ; i ++ )
////    {
////      for( int j = 0 ; j < 2 ; j++ )
////      {
////        BaseComplexNumber.printObject( MultiplyCompMatTest[i][j] ) ;
////      }      
////    }
////    BaseComplexNumber[][] mat1Test = new BaseComplexNumber[2][2] ;
////    ComplexNum temp =  new ComplexNum( 1.0, 1.0 ) ;
////    mat1Test[0][0] = temp ;
////    mat1Test[0][1] = temp ;
////    ComplexNum temp2 =  new ComplexNum( 0.0, 1.0 ) ;
////    mat1Test[1][0] = temp2 ;
////    mat1Test[1][1] = temp2 ;
////    double[] mat2Test = new double[2] ;
////    mat2Test[0] = 2.0 ;
////    mat2Test[1] = 3.0 ;
////    BaseComplexNumber[][] MulColumnMatCompMatTest = new BaseComplexNumber[2][2] ;
////    MulColumnMatCompMatTest = MulColumnMatCompMat( mat2Test, mat1Test  ) ;
////    for( int i = 0 ; i < 2 ; i ++ )
////    {
////      for( int j = 0 ; j < 2 ; j++ )
////      {
////        BaseComplexNumber.printObject( MulColumnMatCompMatTest[i][j] ) ;
////      }      
////    }
////    double[] pa = new double[11] ;
////    pa[0] = 0.0 ;
////    pa[1] = 1.0/100.0 ;
////    pa[2] = 4.0/100.0 ;
////    pa[3] = 9.0/100.0 ;
////    pa[4] = 16.0/100.0 ;
////    pa[5] = 25.0/100.0 ;
////    pa[6] = 36.0/100.0 ;
////    pa[7] = 49.0/100.0 ;
////    pa[8] = 64.0/100.0 ;
////    pa[9] = 81.0/100.0 ;
////    pa[10] = 100.0/100.0 ;
////    double ITest ;
////    ITest = trapezoidalRule(  pa , 10.0 ) ;
////    System.out.println( "ITest" + ITest ) ;
////    int N = 10;
////    double R = 50*Math.pow(10,-3) ;
////    double b = 10*Math.pow(10,-3) ;
//
//    //double a = 3*Math.pow(10,-3) ;
//    //double delta = 0.05*Math.pow(10,-3) ;
//    //double E1 = 2*Math.pow(10,9) ;
//    //double mu1 = 0.3 ;
//
//
//    //BaseComplexNumber[] paResult = new BaseComplexNumber[N+1] ;
//    //paResult = contactDistribution( R , b , P , a , delta , E1 , mu1 , E2 , mu2 , N ) ;
//    //double[] paResulttemp = new double[N+1] ;
//    //for( int i = 0 ; i < N+1 ; i++)
//    //{
//      //paResulttemp[i] = paResult[i].getMod() ; 
//      //System.out.println( paResulttemp[i] ) ;
//    //}
//    
//    double[] surfaceDisp1 = new double[N+1] ;
//    surfaceDisp1 = vpiAlpha( mu1, E1 , b , N , a ) ;
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      System.out.println("i" + i + "surfaceDisp1" + surfaceDisp1[i]) ;
////    }
//    double[] surfaceDisp2 = new double[N+1] ;
//    double[] surfaceDisp2temp = new double[N+1] ;
//    surfaceDisp2temp = vpiAlpha( mu2, E2 , b , N , a ) ;
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      surfaceDisp2[i] = - surfaceDisp2temp[i] ;
//    }
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      System.out.println("i" + i + "surfaceDisp2" + surfaceDisp2[i]) ;
////    }
//    double[][] samplePoints = samplingPoints( N , a ) ;
//    for( int i = 0 ; i < 2 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        System.out.println("i" + i + "j" + j +"samplePoints" + samplePoints[i][j]) ;
//      }  
//    }
//    
//    double[] vy = valpha( mu2, E2 , b , N , a ) ;
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      System.out.println("i" + i + "vy" + vy[i]) ;
//    }
////    BaseComplexNumber[][] fftcoeffN = FastFourierTransform.fftCoeffDomain( samplePoints[0] , samplePoints[1] , N ) ;
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      for( int j = 0 ; j < N+1 ; j++ )
////      {
////        System.out.println("i" + i + "j" + j +"fftcoeffN" + fftcoeffN[i][j]) ;
////      }  
////    }
//    
//    BaseComplexNumber[][] fftcoeffNReal = FastFourierTransform.fftCoeffRealDomain( samplePoints[0] , samplePoints[1] , N ) ;
////  for( int i = 0 ; i < N+1 ; i++ )
////  {
////    for( int j = 0 ; j < N+1 ; j++ )
////    {
////      System.out.println("i" + i + "j" + j +"fftcoeffN" + fftcoeffN[i][j]) ;
////    }  
////  }
////    BaseComplexNumber[][] ifftcoeffN = FastFourierTransform.ifftCoeffDomain( samplePoints[0] , samplePoints[1] , N ) ;
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      for( int j = 0 ; j < N+1 ; j++ )
////      {
////        System.out.println("i" + i + "j" + j +"ifftcoeffN" + ifftcoeffN[i][j]) ;
////      }  
////    }
//    BaseComplexNumber[][] ifftcoeffNReal = FastFourierTransform.ifftCoeffRealDomain( samplePoints[0] , samplePoints[1] , N ) ;
////  for( int i = 0 ; i < N+1 ; i++ )
////  {
////    for( int j = 0 ; j < N+1 ; j++ )
////    {
////      System.out.println("i" + i + "j" + j +"ifftcoeffN" + ifftcoeffN[i][j]) ;
////    }  
////  }
////    BaseComplexNumber[][] T1temp = MulColumnMatCompMat( surfaceDisp1 , fftcoeffN ) ;
//    BaseComplexNumber[][] T1temp = MulColumnMatCompMat( surfaceDisp1 , fftcoeffNReal ) ;
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      for( int j = 0 ; j < N+1 ; j++ )
////      {
////        System.out.println("i" + i + "j" + j +"T1temp" + T1temp[i][j]) ;
////      }  
////    }
////    BaseComplexNumber[][] T1 = MultiplyCompMat( ifftcoeffN , T1temp ) ;
//    BaseComplexNumber[][] T1 = MultiplyCompMat( ifftcoeffNReal , T1temp ) ;
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      for( int j = 0 ; j < N+1 ; j++ )
////      {
////        System.out.println("i" + i + "j" + j +"T1" + T1[i][j]) ;
////      }  
////    }
////    BaseComplexNumber[][] T2temp = MulColumnMatCompMat( surfaceDisp2 , fftcoeffN ) ;
//    BaseComplexNumber[][] T2temp = MulColumnMatCompMat( surfaceDisp2 , fftcoeffNReal ) ;
//    
////    BaseComplexNumber[][] T2 = MultiplyCompMat( ifftcoeffN , T2temp ) ;
//    BaseComplexNumber[][] T2 = MultiplyCompMat( ifftcoeffNReal , T2temp ) ;
//    
//    
////    for( int i = 0 ; i < N+1 ; i++ )
////    {
////      for( int j = 0 ; j < N+1 ; j++ )
////      {
////        System.out.println("i" + i + "j" + j +"T2" + T2[i][j]) ;
////      }  
////    }
//    BaseComplexNumber[][] T = addCompMat( T1 , T2 ) ;
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      for( int j = 0 ; j < N+1 ; j++ )
//      {
//        System.out.println("i" + i + "j" + j +"T" + T[i][j]) ;
//      }  
//    }
//    BaseComplexNumber[] ha = new BaseComplexNumber[N+1] ;
//    ha = gapFunction( R , delta , a , N ) ;
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      //System.out.println("i" + i + "ha" + ha[i]) ;
//    }
//    BaseComplexNumber[] pa = new BaseComplexNumber[N+1] ;
//    pa = GaussElimination.lsolve( T , ha ) ;
//    for( int i = 0 ; i < N+1 ; i++ )
//    {
//      System.out.println("i" + i + "pa" + pa[i]) ;
//    }
//    
//    
////    double[] paMod = new double[pa.length] ;
////    for( int i = 0 ; i < pa.length ; i++ )
////    {
////      paMod[i] = pa[i].getMod() ;
////      System.out.println( "" + paMod[i] ) ;
////    }
////    
//    
    
    // pDistribution( R , b , P , a , delta , E1 , mu1 , E2 , mu2 , N );
    
    // paConverge( R , b , P , a , delta , E1 , mu1 , E2 , mu2 );
    
  }
  

  
    
}

