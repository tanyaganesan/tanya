
public abstract class BaseComplexNumber
{
  public static void printArray( BaseComplexNumber[] baseComplexNumber )
  {
    System.out.println( " printing baseComplexNumber BaseComplexNumber[]" ) ;
    for( int i = 0 ; i < baseComplexNumber.length ; i++ )
    {
      System.out.println( "baseComplexNumber[i]" + baseComplexNumber[i] ) ;
    }
  }
  
  public static void printObject( BaseComplexNumber baseComplexNumber )
  {
    System.out.println( "baseComplexNumber" + baseComplexNumber ) ;
  }
  
  public ComplexNum add( BaseComplexNumber baseComplexNumber )
  {
    double addedReal = baseComplexNumber.getRealPart() + getRealPart() ;
    double addedImaginary = baseComplexNumber.getImaginaryPart() + getImaginaryPart() ;
    ComplexNum added = new ComplexNum( addedReal, addedImaginary ) ;
    return added ;
  }
  public ComplexNum multiply( BaseComplexNumber baseComplexNumber )
  {
//    printObject(  this ) ;
//    printObject(  baseComplexNumber ) ;
//    System.out.println( "baseComplexNumber.getRealPart()" + baseComplexNumber.getRealPart() ) ;
//    System.out.println( "getRealPart()" + getRealPart() ) ;
    double multipliedReal = baseComplexNumber.getRealPart()*getRealPart() - baseComplexNumber.getImaginaryPart()*getImaginaryPart() ;
//    System.out.println( "multipliedReal" + multipliedReal ) ;
    double multipliedImaginary = baseComplexNumber.getRealPart()*getImaginaryPart() + baseComplexNumber.getImaginaryPart()*getRealPart() ;
//    System.out.println( "multipliedImaginary" + multipliedImaginary ) ;
    ComplexNum multiplied = new ComplexNum( multipliedReal, multipliedImaginary ) ;
    return multiplied ; 
  }
  
  public ComplexNum divide( BaseComplexNumber baseComplexNumber )
  {
    ComplexNum conjugate = new ComplexNum( baseComplexNumber.getRealPart(), -baseComplexNumber.getImaginaryPart() ) ;
    ComplexNum nr = multiply(conjugate) ;
    double reciprocalModSqr = 1.0 / Math.pow( baseComplexNumber.getMod() , 2 ) ;
    ComplexNum divided = nr.multiplyConstant(reciprocalModSqr) ;
    return divided ;
  }
  
  public ComplexNum multiplyConstant( double constant )
  {
    double multipliedReal = constant*getRealPart() ;
//    System.out.println( "multipliedReal" + multipliedReal ) ;
    double multipliedImaginary = constant*getImaginaryPart() ;
//    System.out.println( "multipliedImaginary" + multipliedImaginary ) ;
    ComplexNum multiplied = new ComplexNum( multipliedReal, multipliedImaginary ) ;
    return multiplied ;
  }
  public ComplexNum subtract( BaseComplexNumber baseComplexNumber )
  {
    double addedReal = getRealPart() - baseComplexNumber.getRealPart() ;
    double addedImaginary = getImaginaryPart() - baseComplexNumber.getImaginaryPart() ;
    ComplexNum added = new ComplexNum( addedReal, addedImaginary ) ;
    return added ;
  }
  
  public ComplexNum divideInt( int n )
  {
    double divideReal = getRealPart()/n ;
    double divideImaginary = getImaginaryPart()/n ;
    ComplexNum dividedInt = new ComplexNum( divideReal, divideImaginary ) ;
    return dividedInt ;
  }
  

  
  public abstract double getRealPart() ;
  public abstract double getImaginaryPart() ;
  public abstract double getMod() ;
}
