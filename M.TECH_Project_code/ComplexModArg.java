
public class ComplexModArg extends ComplexRootUnity
{
  int r_ ;
  public ComplexModArg( int n , int r )
  {
    super( n ) ;
    r_ = r ;
  }
  
  public double getRealPart()
  {
    return super.getRealPart()*r_ ;
  }
  
  public double getImaginaryPart()
  {
    return super.getImaginaryPart()*r_ ;
  }
  
 
}
