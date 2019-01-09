public class ComplexRootUnity extends BaseComplexNumber 
{
  int n_ = 0 ;
  
  public static final double epsilon_ = 0.0001 ;
  
  public  ComplexRootUnity( int n ) // nth root of unity
  {
    super() ;
    n_ = n ;    
  }
  
  public  double getRealPart()
  {
    double real = Math.cos( 2*Math.PI/n_ ) ;
    if( real < 0.0001 & real > - 0.0001)
    {
      real = 0 ;
    }
    return real ;
  }
  public  double getImaginaryPart()
  {
    double imaginary = Math.sin( 2*Math.PI/n_ ) ;
    if( imaginary < 0.0001 & imaginary > - 0.0001)
    {
      imaginary = 0 ;
    }
    return imaginary  ;
  }
  
  public double getMod()
  {
    double mod ;
    double temp1 = Math.pow(  getRealPart() , 2.0 ) ;
    double temp2 = Math.pow( getImaginaryPart( ) , 2.0 ) ;
    mod = Math.sqrt( temp1 + temp2 ) ;
    return mod ;
  }
  
  
  public String toString()
  {
    return ( super.toString() + " n_ " + n_ + " realRootUnity" + getRealPart() + " imaginaryRootUnity" + getImaginaryPart() ) ; 
  }
 
  
}
