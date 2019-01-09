
public class ComplexRootArg extends BaseComplexNumber
{
  double arg_ ;
  
  public static final double epsilon_ = 0.0001 ;
  
  public  ComplexRootArg( double arg ) // nth root of unity
  {
    super() ;
    arg_ = arg ;    
  }
  
  public  double getRealPart()
  {
    double real = Math.cos( arg_ ) ;
    if( real < 0.0001 & real > - 0.0001)
    {
      real = 0 ;
    }
    return real ;
  }
  
  public  double getImaginaryPart( )
  {
    double imaginary = Math.sin( arg_ ) ;
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
    return ( super.toString() + " real " + getRealPart() + " imaginary " + getImaginaryPart( ) ) ; 
  }
}
