
public class ComplexNum extends BaseComplexNumber
{
  double real_ ;
  double imaginary_ ;
  
  public ComplexNum( double real, double imaginary )
  {
    super() ;
    real_ = real ;
    imaginary_ = imaginary ;
  }
  
  public double getRealPart()
  {
    return real_ ;
  }

  public double getImaginaryPart()
  {
    return imaginary_ ;
  }
  
  
  public double getMod()
  {
    double mod ;
    double temp1 = Math.pow( real_ , 2.0 ) ;
    double temp2 = Math.pow( imaginary_ , 2.0 ) ;
    mod = Math.sqrt( temp1 + temp2 ) ;
    return mod ;
  }
  
  public String toString()
  {
    return ( super.toString() + " real_ " + real_ + " imaginary_  " + imaginary_ ) ; 
  }
}
