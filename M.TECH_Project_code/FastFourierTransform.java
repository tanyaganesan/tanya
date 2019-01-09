
public class FastFourierTransform
{
  /////Given signal is periodic with 2a. Bandwidth 2B. N sampling points. Domain 0 to 2a. Frequency domain 0 to 2B
  public static BaseComplexNumber[][] fftCoeff( int N )
  {
    BaseComplexNumber[][] fftCoeff = new BaseComplexNumber[N][N] ;
    for( int i = 0 ; i < N ; i++ )
    {
      for( int j = 0 ; j < N ; j++ )
      {
        double argTemp = 2.0* Math.PI * i * j / N ;
        ComplexRootArg temp = new ComplexRootArg( argTemp ) ;
        fftCoeff[i][j] = temp ;
      }
    }
    return fftCoeff ;
  }
  
  public static BaseComplexNumber[][] fftCoeffDomain( double[] x , double[] alpha , int N)
  {
    ////// x domain , frequency domain and number of sampling points given.
    BaseComplexNumber[][] fftCoeff = new BaseComplexNumber[ N + 1 ][ N + 1 ] ;
    for( int i = 0 ; i < N + 1 ; i++ )
    {
      for( int j = 0 ; j < N + 1 ; j++ )
      {
        double argTemp =  alpha[i]*x[j] ;
        ComplexRootArg temp = new ComplexRootArg( argTemp ) ;
        fftCoeff[i][j] = temp ;
      }
    }
    return fftCoeff ;
  }
  
  ////////Fourier cosine terms. Used in case the function is real.
  public static BaseComplexNumber[][] fftCoeffRealDomain( double[] x , double[] alpha , int N)
  {
    ////// x domain , frequency domain and number of sampling points given.
    BaseComplexNumber[][] fftCoeffReal = new BaseComplexNumber[ N + 1 ][ N + 1 ] ;
    for( int i = 0 ; i < N + 1 ; i++ )
    {
      for( int j = 0 ; j < N + 1 ; j++ )
      {
        double argTemp =  alpha[i]*x[j] ;
        fftCoeffReal[i][j] = new ComplexNum( Math.cos( argTemp ) , 0.0 ) ;
      }
    }
    return fftCoeffReal ;
  }
  
  public static BaseComplexNumber[][] ifftCoeff( int N )
  {
    BaseComplexNumber[][] ifftCoeff = new BaseComplexNumber[N][N] ;
    for( int i = 0 ; i < N ; i++ )
    {
      for( int j = 0 ; j < N ; j++ )
      {
        double argTemp = -2.0* Math.PI * i * j / N ;
        ComplexRootArg temp = new ComplexRootArg( argTemp ) ;
        ifftCoeff[i][j] = temp.divideInt( N ) ;
      }
    }
    return ifftCoeff ;
  }
  
  public static BaseComplexNumber[][] ifftCoeffDomain( double[] x , double[] alpha , int N)
  {
    ////// x domain , frequency domain and number of sampling points given.
    BaseComplexNumber[][] ifftCoeff = new BaseComplexNumber[N+1][N+1] ;
    for( int i = 0 ; i < N+1 ; i++ )
    {
      for( int j = 0 ; j < N+1 ; j++ )
      {
        double argTemp =  -alpha[i]*x[j] ;
        ComplexRootArg temp = new ComplexRootArg( argTemp ) ;
        ifftCoeff[i][j] = temp ;
      }
    }
    return ifftCoeff ;
  }
  
  ////////Fourier inverse cosine terms. Used in case the function is real.
  public static BaseComplexNumber[][] ifftCoeffRealDomain( double[] x , double[] alpha , int N)
  {
    ////// x domain , frequency domain and number of sampling points given.
    BaseComplexNumber[][] ifftCoeffReal = new BaseComplexNumber[ N + 1 ][ N + 1 ] ;
    for( int i = 0 ; i < N + 1 ; i++ )
    {
      for( int j = 0 ; j < N + 1 ; j++ )
      {
        double argTemp =  -alpha[i]*x[j] ;
        ifftCoeffReal[i][j] = new ComplexNum( Math.cos(argTemp), 0.0 ) ;
      }
    }
    return ifftCoeffReal ;
  }
  
  //N can be only exponential of 2 
  public static BaseComplexNumber[] fft( BaseComplexNumber[] ba )
  {
    int n = ba.length ;
    BaseComplexNumber[] y = new ComplexNum[n] ;
    if( n == 1 )
    {
      y = ba ;
      return y ;
    }
    else
    {
      BaseComplexNumber omegaN = new ComplexRootUnity( n ) ;
//      System.out.println("OmegaN");
//      BaseComplexNumber.printObject( omegaN ) ;
      BaseComplexNumber omega = new ComplexRootUnity( 1 ) ;
//    System.out.println("Omega");
//    BaseComplexNumber.printObject( omega ) ;
      BaseComplexNumber[] ba0 = new BaseComplexNumber[n/2] ;
      BaseComplexNumber[] ba1 = new BaseComplexNumber[n/2] ;
      
      int temp1 = 0 ;
      for( int i = 0 ; i < n/2 ; i++ )
      {
        ba0[i] = ba[temp1] ;
        temp1 = temp1 + 2 ;
      }
      
      int temp2 = 1 ;
      for( int j = 0 ; j < n/2 ; j++ )
      {
        ba1[j] = ba[temp2] ;
        temp2 = temp2 + 2 ;
      }
      
      BaseComplexNumber[] y0 = new BaseComplexNumber[n/2] ;
      y0 = fft( ba0 ) ;
//      System.out.println( "y0") ;
//      BaseComplexNumber.printArray( y0 );
      
      BaseComplexNumber[] y1 = new BaseComplexNumber[n/2] ;
      y1 = fft( ba1 ) ;
//      System.out.println( "y1") ;
//      BaseComplexNumber.printArray( y1 );
      
      for( int k = 0 ; k <= n/2 - 1 ; k++ )
      {
        BaseComplexNumber temp3 = omega.multiply( y1[k] ) ;
//        System.out.println("Here0");
//        BaseComplexNumber.printObject( temp3 ) ;
        y[k] = y0[k].add( temp3 ) ;
//        System.out.println("Here1");
//        BaseComplexNumber.printObject( y[k] ) ;
        y[k + n/2] = y0[k].subtract( temp3 ) ;
//        BaseComplexNumber.printObject( y[k + n/2] ) ;
//        System.out.println("Here2");
        omega = omega.multiply( omegaN ) ;
      }
    }
    return y ;
  }
  
  public static BaseComplexNumber[] ifft( BaseComplexNumber[] by )
  {
    int n = by.length ;
    BaseComplexNumber[] ba = new ComplexNum[n] ;
    if( n == 1 )
    {
      ba = by ;
      return ba ;
    }
    else
    {
      BaseComplexNumber omegaN = new ComplexRootUnity( - n ) ;
      BaseComplexNumber omega = new ComplexRootUnity( 1 ) ;
      BaseComplexNumber[] by0 = new BaseComplexNumber[n/2] ;
      BaseComplexNumber[] by1 = new BaseComplexNumber[n/2] ;
      
      int temp1 = 0 ;
      for( int i = 0 ; i < n/2 ; i++ )
      {
        by0[i] = by[temp1] ;
        temp1 = temp1 + 2 ;
      }
      
      int temp2 = 1 ;
      for( int j = 0 ; j < n/2 ; j++ )
      {
        by1[j] = by[temp2] ;
        temp2 = temp2 + 2 ;
      }
      
      BaseComplexNumber[] a0 = new BaseComplexNumber[n/2] ;
      a0 = ifft( by0 ) ;
      
      BaseComplexNumber[] a1 = new BaseComplexNumber[n/2] ;
      a1 = ifft( by1 ) ;
      
      for( int k = 0 ; k <= n/2 - 1 ; k++ )
      {
        BaseComplexNumber temp3 = omega.multiply( a1[k] ) ;
        ba[k] = a0[k].add( temp3 ) ;
        ba[k + n/2] = a0[k].subtract( temp3 ) ;
        omega = omega.multiply( omegaN ) ;
      }
    }
    return ba ;
  }
    
  public static void lowPass( double[] fs, double[] sa ,double mu0 )
  {
    int N = fs.length ;
    double[] ha = new double[N] ; 
    double[] fh = new double[N] ;
    for( int i = 0 ; i < N ; i++)
    {
      if( sa[i]<mu0 && sa[i]>-mu0 )
      {
        ha[i] = 1 ;
      }
      else
      {
        ha[i] = 0 ;
      }
      fh[i] = ha[i]*fs[i] ;
    }
  }
    public static void main( String[] arg )
    {
//      int n = 4 ;
//      BaseComplexNumber[] ba = new ComplexNum[n] ;
//      ComplexNum ba0 = new ComplexNum( 0.0, 0.0 ) ;
//      ComplexNum ba1 = new ComplexNum( 1.0, 0.0 ) ;
//      ComplexNum ba2 = new ComplexNum( 1.0, 0.0 ) ;
//      ComplexNum ba3 = new ComplexNum( 0.0, 0.0 ) ;
//      ba[0] = ba0 ;
//      ba[1] = ba1 ; 
//      ba[2] = ba2 ;
//      ba[3] = ba3 ;
//      
//      BaseComplexNumber[] result1 = new BaseComplexNumber[n] ;
//      result1 = fft(ba) ;
//      
//      BaseComplexNumber.printArray( result1 );
//      
//      BaseComplexNumber[] result2 = new BaseComplexNumber[n] ;
//      result2 = ifft( result1 ) ;
//      for( int l = 0 ; l < n ; l++ )
//      {
//        result2[l] = result2[l].divideInt( n ) ;
//      }
//      BaseComplexNumber.printArray( result2 ); 
      
//      int N = 4 ;
//      BaseComplexNumber[][] fftCoeffResult = new BaseComplexNumber[N][N] ;
//      fftCoeffResult = fftCoeff( N ) ;
//      for( int i = 0 ; i < N ; i++ )
//      {
//        for( int j = 0 ; j < N ; j++ )
//        {
//          System.out.println( "fftCoeffResult[i][j]" + fftCoeffResult[i][j] ) ;
//        }
//      }
//      
//    int N = 4 ;
//    BaseComplexNumber[][] ifftCoeffResult = new BaseComplexNumber[N][N] ;
//    ifftCoeffResult = ifftCoeff( N ) ;
//    for( int i = 0 ; i < N ; i++ )
//    {
//      for( int j = 0 ; j < N ; j++ )
//      {
//        System.out.println( "fftCoeffResult[i][j]" + ifftCoeffResult[i][j] ) ;
//      }
//    }
//      
      /////////// Fourier transform of sine wave
//      int N = 512 ;
//      BaseComplexNumber[] cosxa = new BaseComplexNumber[N] ;
//      double[] xa = new double[N] ;
//      xa[0] = 0.0 ;
//      for( int i = 0 ; i < N ; i++ )
//      {
//        BaseComplexNumber y = new ComplexNum( Math.cos( 2.0*2.0*Math.PI*xa[i] ) , 0.0 ) ;
//        cosxa[i] = y ;
//        if( i < N -1 )
//        {
//          xa[i+1] = xa[i] + 0.1 ;
//        }        
//      }
//      double bandwidth = N/0.5 ;
//      System.out.println(bandwidth) ;
      
//      for( int i = 0 ; i < N ; i++)
//      {
//        System.out.println(sinxa[i]) ;
//      }
      
//      BaseComplexNumber[] fs = new BaseComplexNumber[N] ;
//      fs = fft(cosxa) ;
//    for( int i = 0 ; i < N ; i++)
//    {
//      System.out.println(fs[i]) ;
//    }
//    
//      
    }
  }



