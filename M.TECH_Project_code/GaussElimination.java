
public class GaussElimination 
{
  public static BaseComplexNumber[] lsolve(BaseComplexNumber[][] A, BaseComplexNumber[] b)
  {
    int N = b.length;
    for (int p = 0; p < N; p++)
    {
      // find pivot row and swap
      int max = p;
      for (int i = p + 1; i < N; i++)
      {
        if ((A[i][p].getMod()) > (A[max][p].getMod()))
        {
          max = i;
        }
      }
      BaseComplexNumber[] temp = A[p];
      A[p] = A[max];
      A[max] = temp;
      BaseComplexNumber t = b[p];
      b[p] = b[max];
      b[max] = t;
      // pivot within A and b
      for (int i = p + 1; i < N; i++)
      {       
        BaseComplexNumber alpha = A[i][p].divide( A[p][p] )  ;
        BaseComplexNumber temp1 = alpha.multiply( b[p] ) ;
        b[i] = b[i].subtract( temp1 )  ;       
        for (int j = p; j < N; j++)
        {
          BaseComplexNumber temp2 = alpha.multiply( A[p][j] ) ;
          A[i][j] = A[i][j].subtract( temp2 )  ;
        }
      }
    }  
    // back substitution
    BaseComplexNumber[] x = new BaseComplexNumber[N];
    for (int i = N - 1; i >= 0; i--)
    {
      BaseComplexNumber sum = new ComplexNum( 0.0, 0.0);
      for (int j = i + 1; j < N; j++)
      {
        BaseComplexNumber temp3 = A[i][j].multiply( x[j] ) ;
        sum = sum.add( temp3 )  ;
      }
      BaseComplexNumber temp4 = b[i].subtract( sum ) ;    
      x[i] = temp4.divide(  A[i][i] );
    }
    
    return x ;   
    }
  
  
  public static void main( String[] arg )
  {
    BaseComplexNumber[][] GaussTestA = new BaseComplexNumber[2][2] ;
    BaseComplexNumber[] GaussTestb = new BaseComplexNumber[2] ;
    BaseComplexNumber[]  GaussTestResult = new BaseComplexNumber[2] ;
    GaussTestA[0][0] = new ComplexNum(1.0, 0.0) ;
    GaussTestA[0][1] = new ComplexNum(1.0, 0.0) ;
    GaussTestA[1][0] = new ComplexNum(-1.0, 1.0) ;
    GaussTestA[1][1] = new ComplexNum(1.0, 0.0) ;
    GaussTestb[0] = new ComplexNum(4.0, 0.0) ;
    GaussTestb[1] = new ComplexNum(-1.0, 0.0) ;
    
    GaussTestResult = lsolve( GaussTestA ,GaussTestb ) ;
    BaseComplexNumber.printArray(GaussTestResult) ;
  }
  
  }
