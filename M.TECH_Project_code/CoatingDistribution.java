
public class CoatingDistribution 
{
	public static double simpsonRule( double x0, double x2 )
	{
	  /////////////// Only for standard normal distribution
	  double h = 0.001 ;
	  double n = ( x2 - x0 )/h ; 
	  //System.out.println("n" + n) ;
	  int size = (int)n ;
	  //System.out.println("size" + size ) ;
	  double[] xa = new double[size] ;
	  double Ix0x2 = 0.0 ;
	  xa[0] = x0 ;
	  for( int i = 1 ; i < n - 1 ; i++ )
	  {
	    xa[i] = xa[i-1] + h ;  
	  }
	  for( int i = 0 ; i < n - 2 ; i++ )
	  {
		 double Ixi = ( 1.0/3.0 )*h*( pdf( xa[i] ) + 4*pdf( ( xa[i] + xa[i+1] )/2  ) + pdf( xa[i+1] )  ) ;
		 //System.out.println("Ixi" + Ixi ) ;
		 Ix0x2 = Ix0x2 + Ixi ;
	  }	  
	  return Ix0x2 ;
	}
	
	public static double pdf( double x )
	{
		double mu = 0.0 ; 
		double sigma = 1.0 ;
		double pdf = ( ( 1.0 / ( Math.pow( 2* Math.PI*Math.pow( sigma, 2.0 ) , 0.5 ) ) )*Math.exp( - Math.pow(( ( x - mu) /sigma ), 2.0) ) ) ;
		//System.out.println("pdf" + pdf ) ;
		return pdf ;
	}
	public static double[][] inclusionsDistribution( double a , double b )
	{
	  //double a = 4.287*Math.pow(10.0,-7.0) ;
	  //double b = 1.90057*Math.pow(10.0,-5.0);
	  
	  int nSpheres = 100 ;
	  int interval = 20 ;
	  double window = (b - a)/(interval) ;
	  double[] ra = new double[interval+1] ;
	  double[][] result = new double[2][interval] ;
	  //////////radii array
	  ra[0] = a ;
	  for( int i = 1 ; i < interval + 1 ; i++ )
	  {
	    ra[i] = ra[i-1] + window ;
	  }
	  //////////mean radii array and volume 
	  double[] raMean = new double[interval] ;
	  double[] volumeSpherea = new double[interval] ;
	  for( int k = 0 ; k < interval ; k++ )
	  {
	    raMean[k] = ( ( ra[k] + ra[k+1] )/2 ) / Math.pow( 10, -6.0 ) ; // in micrometer
	    volumeSpherea[k] = (4.0/3.0)*Math.PI*Math.pow( raMean[k], 3 ) ;
	  }
	  result[0] = raMean ;
	  ///////////mean 
	  double sum_r = 0 ;
	  for( int j = 0 ; j < interval + 1 ; j++ )
	  {
	    sum_r = sum_r + ra[j] ;
	  }
	  double mu = sum_r / ( interval + 1 ) ;
	  ///////standard deviation
	  double sigma_sum_r = 0 ;
	  for( int j = 0 ; j < interval + 1 ; j++ )
	  {
		  sigma_sum_r = sigma_sum_r + Math.pow( ra[j] - mu , 2 ) ;
	  }
	  double sigma = Math.sqrt(sigma_sum_r / (interval + 1)) ;
	  ///////probability array
	  double z1 ;
	  double z2 ;
	  double[] pa = new double[interval] ;
	  double[] spherea = new double[interval] ;
	  //double temp = 0.0 ;
	  for( int l = 0 ; l < interval ; l++ )
	  {
	    z1 = ( ra[l] - mu ) / sigma ; 
	    //System.out.println("z1" + z1) ;
	    z2 = ( ra[l+1] - mu ) / sigma ;
	    //System.out.println("z2" + z2) ;
	    pa[l] = simpsonRule( z1 , z2 ) ;
	    //System.out.println("pa[l]" + pa[l]) ;
	    //temp = temp + pa[l] ;
	    spherea[l] = pa[l] * nSpheres ;
	  }
	  //System.out.println("temp" + temp) ;
	  //////volume of spheres and sum of volume 
	  double[] TotalVolumeSpherea = new double[interval] ;
	  double sumVolume = 0.0 ;
	  for( int m = 0 ; m < interval ; m++ )
	  {
		  TotalVolumeSpherea[m] = spherea[m]*volumeSpherea[m] ;
		  sumVolume = sumVolume + TotalVolumeSpherea[m] ; 	 
	  }
	  //System.out.println( "sumVolume" + sumVolume ) ;
	  //////Partial Volume fractions
	  double[] pvfSpherea= new double[interval] ;
	  double vf = 0.1 ;
	  for( int n = 0 ; n < interval ; n++ )
	  {
		  pvfSpherea[n] = ( TotalVolumeSpherea[n] / sumVolume )*100.0*vf ;
		  //System.out.println( "pvfSpherea" + pvfSpherea[n] ) ;
	  }	 
	  result[1] = pvfSpherea ;
	  return result ;
 	}
	

	public static void main( String[] arg )
	{
	  double a = 0.4287*Math.pow(10.0,-7.0) ;
	  double b = 0.190057*Math.pow(10.0,-5.0);

//	  for( int i = 0 ; i < pvf.length ; i++ )
//	  {
//        System.out.println( "pvf" + pvf[i] ) ;
//	  }
	  double side = (1.429*Math.pow(10.0,-5.0))/ Math.pow( 10, -6.0 ) ; // in micrometer 
	  double vf_tolerance = 0.0001 ;
	  double n = 999.0 ;// constant
	  int interval = 20 ;
	  StringBuffer sb = new StringBuffer() ;
	  sb.append("thesisoutd") ;
	  sb.append('\n') ;
	  sb.append(0) ;
	  sb.append('\t') ;
	  sb.append("100000000") ;
	  sb.append('\t') ;
	  sb.append("1") ;
	  sb.append('\t') ;
	  sb.append(interval) ;
	  sb.append('\t') ;
	  sb.append('\t') ;
	  sb.append("201") ;
	  sb.append('\n') ;
	  sb.append(side) ;
	  sb.append(n) ;
	  sb.append('\t') ;
	  sb.append(vf_tolerance) ;
	  sb.append('\n') ;
	  double[][] meanRpvf = inclusionsDistribution(a,b) ;
	  for( int i = 0 ; i < meanRpvf.length ; i++)
	  {
	    for( int j = meanRpvf[1].length-1 ; j >= 0 ; j-- )
	    {
	    	sb.append( meanRpvf[i][j] ) ;
	    	sb.append(' ') ;
	    }
	    sb.append('\n') ;
	  }
//	  double window = (b - a)/(interval) ;
//	  for( double i = 0.0 ; i < interval + 1 ; i++ )
//	  {
//	    double temp = a + i*window ;
//	    sb.append(temp) ;
//	    sb.append(" ") ;
//	  }
//	  double mu = 0.0 ;
//	  double sigma = 1.0 ;
//	  double x ;
//	  double pdf = ( ( 1 / ( Math.pow( 2* Math.PI*Math.pow( sigma, 2.0 ) , 0.5 ) ) )*Math.exp( - Math.pow(( ( x - mu) /sigma ), 2.0) ) ) ;
//	  double x0 ;
//	  double x1 ;
//	  double x2 ;
//	  double Ix1x2 ;
//	  double h ;
//	  Ix1x2 = 
	  System.out.print(sb.toString());
	}
	
}

