package Entities;

import java.util.ArrayList;
import java.util.List;

import org.lwjgl.util.vector.Vector3f;
import org.lwjgl.util.vector.Vector4f;


public class Universe {
	private List<Unit> units = new ArrayList<Unit>();
	private final static float GRAVITATIONAL_CONSTANT = (float) (11.674 * Math.pow(10, -9));

	public Universe(){		
	}

	public void updateUniverse(Player player){
		for (int i = 0; i < units.size(); i++) {
			units.get(i).update();
		}
		player.update();
		detectCollision(player);
	}

	public void addUnit(Unit unit){
		unit.setUnitID(units.size());
		this.units.add(unit);
	}

	public void removeUnit(int unitID){
		this.units.remove(unitID);
	}

	public List<Unit> getUniverse(){
		return units;
	}

	public Vector3f findPointOnLine(Vector3f point1, Vector3f point2, float dist){
		float dx, dy, dz;

		dx = point2.x - point1.x;
		dy = point2.y - point1.y;
		dz = point2.z - point1.z;
		
		return new Vector3f((dx * dist), (dy * dist), (dz * dist));
	}

	public void detectCollision(Player player){
		//GL20.glGetVertexAttribPointer(0, GL20., 4);
		int x;
		Unit currentUnit;
		Unit tempUnit;
		for (x = 0; x < units.size()+1; x++) {
			if (x == units.size()){
				currentUnit = player;
			} else {
				currentUnit = units.get(x);
			}
			for (int y = 0; y < units.size()+1; y++) {
				if (y != x){
					if (y == units.size()){
						tempUnit = player;
					} else {
						tempUnit = units.get(y);
					}
					Vector3f offset = findDifference(currentUnit.getPosition(), tempUnit.getPosition());
					float distance = (float) Math.sqrt(Math.pow(offset.x, 2)+Math.pow(offset.y, 2)+Math.pow(offset.z, 2));
					if (distance > (currentUnit.getScale()+tempUnit.getScale())){
						float accelleration = (float) (((GRAVITATIONAL_CONSTANT*((tempUnit.getMass()*currentUnit.getMass())/Math.pow(distance, 2)))/tempUnit.getMass()));
						Vector3f changeVector = findPointOnLine(currentUnit.getPosition(),tempUnit.getPosition(), accelleration);
						currentUnit.increaseDirectionalSpeed(changeVector.x, changeVector.y, changeVector.z);
					} else {
						int error = 0;
						error = collision3D(0.5f, currentUnit, tempUnit, new Vector3f(tempUnit.getPosition().x+(offset.x/2),tempUnit.getPosition().y+(offset.y/2),tempUnit.getPosition().z+(offset.z/2)), error);
						if (error !=0){
							//System.out.println(error);
						}
					}
				}
			}
		}

	}

	public void compareUnits(Unit unit1, Unit unit2){
		Vector3f[] unit1Points = unit1.getUnitPoints();
		Vector3f[] unit2Points = unit2.getUnitPoints();
		List<Vector3f> intersections = new ArrayList<Vector3f>();
		int[] pattern1 = {0,0,0,1};
		int[] pattern2 = {1,2,1,2};
		int[] pattern3 = {2,3,3,3};
		for (int i = 0; i < 4; i++){
			if (pattern1[i] != 0 && pattern2[i] != 1 && pattern3[i] != 2){
				intersections.add(comparePlanes(unit1Points[0],unit1Points[1], unit1Points[2],unit2Points[pattern1[i]],unit2Points[pattern2[i]], unit2Points[pattern3[i]]));
			} 
			if (pattern1[i] != 0 && pattern2[i] != 2 && pattern3[i] != 3){
				intersections.add(comparePlanes(unit1Points[0],unit1Points[2], unit1Points[3],unit2Points[pattern1[i]],unit2Points[pattern2[i]], unit2Points[pattern3[i]]));
			}
			if (pattern1[i] != 0 && pattern2[i] != 1 && pattern3[i] != 3){
				intersections.add(comparePlanes(unit1Points[0],unit1Points[1], unit1Points[3],unit2Points[pattern1[i]],unit2Points[pattern2[i]], unit2Points[pattern3[i]]));
			} 
			if (pattern1[i] != 1 && pattern2[i] != 2 && pattern3[i] != 3){
				intersections.add(comparePlanes(unit1Points[1],unit1Points[2], unit1Points[3],unit2Points[pattern1[i]],unit2Points[pattern2[i]], unit2Points[pattern3[i]]));
			}
		}
		for (int i = 0; i < intersections.size(); i++){
			if (intersections.get(i) != null){
				int error = 0;
				error = collision3D(1f, unit1, unit2, intersections.get(i), error);
				if (error != 0){
					System.out.println("error "+error+" has occurred");
				}
			}
		}
	}

	//  The Parameters are:
	//
	//		    R    (restitution coefficient)  between 0 and 1 (1=perfectly elastic collision)
	//		    m1    (mass of ball 1)
	//		    m2    (mass of ball 2)
	//		    r1    (radius of ball 1)
	//		    r2    (radius of ball 2)
	//  & x1    (x-coordinate of ball 1) 
	//  & y1    (y-coordinate of ball 1)          
	//  & z1    (z-coordinate of ball 1) 
	//  & x2    (x-coordinate of ball 2)              
	//  & y2    (y-coordinate of ball 2)         
	//  & z2    (z-coordinate of ball 2)         
	//  & vx1   (velocity x-component of ball 1) 
	//  & vy1   (velocity y-component of ball 1)
	//  & vz1   (velocity z-component of ball 1)          
	//  & vx2   (velocity x-component of ball 2)         
	//  & vy2   (velocity y-component of ball 2)
	//  & vz2   (velocity z-component of ball 2)
	//  & error (int)     (0: no error 
	//		                     1: balls do not collide
	//		                     2: initial positions impossible (balls overlap))

	public int collision3D(float R, Unit unit1, Unit unit2, Vector3f collisionPoint, int error){ 

		float  pi,r12,m21,d,v,theta2,phi2,st,ct,sp,cp,vx1r,vy1r,vz1r,fvz1r,
		thetav,phiv,dr,alpha,beta,sbeta,cbeta,t,a,dvz2,
		vx2r,vy2r,vz2r,x21,y21,z21,vx21,vy21,vz21,vx_cm,vy_cm,vz_cm;

		float m1 = unit1.getMass();
		float m2 = unit2.getMass();
		Vector3f r1d = findDifference(unit1.getPosition(), collisionPoint);
		float r1 = (float) Math.sqrt(Math.pow(r1d.x, 2)+Math.pow(r1d.y, 2)+Math.pow(r1d.z, 2));
		Vector3f r2d = findDifference(unit2.getPosition(), collisionPoint);
		float r2 = (float) Math.sqrt(Math.pow(r2d.x, 2)+Math.pow(r2d.y, 2)+Math.pow(r2d.z, 2));
		float x1 = unit1.getPosition().x;
		float y1 = unit1.getPosition().y;
		float z1 = unit1.getPosition().z;
		float x2 = unit2.getPosition().x;
		float y2 = unit2.getPosition().y;
		float z2 = unit2.getPosition().z;
		float vx1 = unit1.getDirectionalSpeed().x;
		float vy1 = unit1.getDirectionalSpeed().y;
		float vz1 = unit1.getDirectionalSpeed().z;
		float vx2 = unit2.getDirectionalSpeed().x; 
		float vy2 = unit2.getDirectionalSpeed().y; 
		float vz2 = unit2.getDirectionalSpeed().z;

		//	     **** initialize some variables ****

		pi=(float) Math.acos(-1.0E0);
		error=0;
		r12=r1+r2;
		m21=m2/m1;
		x21=x2-x1;
		y21=y2-y1;
		z21=z2-z1;
		vx21=vx2-vx1;
		vy21=vy2-vy1;
		vz21=vz2-vz1;

		vx_cm = (m1*vx1+m2*vx2)/(m1+m2) ;
		vy_cm = (m1*vy1+m2*vy2)/(m1+m2) ;
		vz_cm = (m1*vz1+m2*vz2)/(m1+m2) ;  


		//	     **** calculate relative distance and relative speed ***
		d=(float) Math.sqrt(x21*x21 +y21*y21 +z21*z21);
		v=(float) Math.sqrt(vx21*vx21 +vy21*vy21 +vz21*vz21);

		//	     **** return if distance between balls smaller than sum of radii ****
		if (d<r12) {error=2; return error;}

		//	     **** return if relative speed = 0 ****
		if (v==0) {error=1; return error;}


		//	     **** shift coordinate system so that ball 1 is at the origin ***
		x2=x21;
		y2=y21;
		z2=z21;

		//	     **** boost coordinate system so that ball 2 is resting ***
		vx1=-vx21;
		vy1=-vy21;
		vz1=-vz21;

		//	     **** find the polar coordinates of the location of ball 2 ***
		theta2=(float) Math.acos(z2/d);
		if (x2==0 && y2==0) {phi2=0;} else {phi2=(float) Math.atan2(y2,x2);}
		st=(float) Math.sin(theta2);
		ct=(float) Math.cos(theta2);
		sp=(float) Math.sin(phi2);
		cp=(float) Math.cos(phi2);


		//	     **** express the velocity vector of ball 1 in a rotated coordinate
		//	          system where ball 2 lies on the z-axis ******
		vx1r=ct*cp*vx1+ct*sp*vy1-st*vz1;
		vy1r=cp*vy1-sp*vx1;
		vz1r=st*cp*vx1+st*sp*vy1+ct*vz1;
		fvz1r = vz1r/v ;
		if (fvz1r>1) {fvz1r=1;}   // fix for possible rounding errors
		else if (fvz1r<-1) {fvz1r=-1;} 
		thetav=(float) Math.acos(fvz1r);
		if (vx1r==0 && vy1r==0) {phiv=0;} else {phiv=(float) Math.atan2(vy1r,vx1r);}


		//	     **** calculate the normalized impact parameter ***
		dr=(float) (d*Math.sin(thetav)/r12);


		//	     **** return old positions and velocities if balls do not collide ***
		if (thetav>pi/2 || Math.abs(dr)>1) {
			x2=x2+x1;
			y2=y2+y1;
			z2=z2+z1;
			vx1=vx1+vx2;
			vy1=vy1+vy2;
			vz1=vz1+vz2;
			error=1;
			return error;
		}

		//	     **** calculate impact angles if balls do collide ***
		alpha=(float) Math.asin(-dr);
		beta=phiv;
		sbeta=(float) Math.sin(beta);
		cbeta=(float) Math.cos(beta);


		//	     **** calculate time to collision ***
		t=(float) ((d*Math.cos(thetav) -r12*Math.sqrt(1-dr*dr))/v);


		//	     **** update positions and reverse the coordinate shift ***
		x2=x2+vx2*t +x1;
		y2=y2+vy2*t +y1;
		z2=z2+vz2*t +z1;
		x1=(vx1+vx2)*t +x1;
		y1=(vy1+vy2)*t +y1;
		z1=(vz1+vz2)*t +z1;



		//  ***  update velocities ***

		a=(float) Math.tan(thetav+alpha);

		dvz2=2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21));

		vz2r=dvz2;
		vx2r=a*cbeta*dvz2;
		vy2r=a*sbeta*dvz2;
		vz1r=vz1r-m21*vz2r;
		vx1r=vx1r-m21*vx2r;
		vy1r=vy1r-m21*vy2r;


		//	     **** rotate the velocity vectors back and add the initial velocity
		//	           vector of ball 2 to retrieve the original coordinate system ****

		vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
		vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
		vz1=ct*vz1r-st*vx1r               +vz2;
		vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
		vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
		vz2=ct*vz2r-st*vx2r               +vz2;


		//	     ***  velocity correction for inelastic collisions ***

		vx1=(vx1-vx_cm)*R + vx_cm;
		vy1=(vy1-vy_cm)*R + vy_cm;
		vz1=(vz1-vz_cm)*R + vz_cm;
		vx2=(vx2-vx_cm)*R + vx_cm;
		vy2=(vy2-vy_cm)*R + vy_cm;
		vz2=(vz2-vz_cm)*R + vz_cm;  

		unit1.setDirectionalSpeed(new Vector3f(vx1, vy1, vz1));
		unit2.setDirectionalSpeed(new Vector3f(vx2, vy2, vz2));
		return error;
	}

	public Vector3f findDifference(Vector3f point1, Vector3f point2){
		Vector3f difference = new Vector3f(0,0,0);
		difference.x = point1.x - point2.x;
		difference.y = point1.y - point2.y;
		difference.z = point1.z -point2.z;
		return difference;
	}

	private Vector3f crossProduct(Vector3f a, Vector3f b){
		//System.out.println(a);
		//System.out.println(b);
		float i = ((a.y*b.z)-(a.z*b.y));
		float j = -((a.x*b.z)-(a.z*b.x));
		float k = ((a.y*b.x)-(a.x*b.y));
		//System.out.println(new Vector3f(i,j,k));
		return new Vector3f(i,j,k);
	}

	private Vector4f formatPlane(Vector3f plane1, Vector3f crossProduct){
		//System.out.println(plane1);
		float constant = 0;
		float i = crossProduct.x;
		constant += plane1.x*crossProduct.x;
		float j = crossProduct.y;
		constant += plane1.y*crossProduct.y;
		float k = crossProduct.z;
		constant += plane1.z*crossProduct.z;
		return new Vector4f(i,j, k, constant);
	}

	private Vector3f intersection(Vector4f planeA, Vector4f planeB){
		//System.out.println(planeA);
		//System.out.println(planeB);
		Vector4f tempPlaneA = new Vector4f(0,0,0,0);
		Vector4f tempPlaneB = new Vector4f(0,0,0,0);
		Vector3f tempPlaneC = new Vector3f(0,0,0);
		Vector3f tempPlaneD = new Vector3f(0,0,0);
		Vector3f tempPlaneX = new Vector3f(0,0,0);
		float j = 0;
		float k = 0;
		float l = 0;
		float x = 0;
		float y = 0;
		float z = 0;
		float factor = planeA.x*planeB.x;
		float factorA = factor/planeA.x;
		float factorB = factor/planeB.x;
		if (factor == 0){
			factorA = 1;
			factorB = 1;
		}
		tempPlaneA = new Vector4f(planeA.x*factorA, planeA.y*factorA, planeA.z*factorA, planeA.w*factorA);
		tempPlaneB = new Vector4f(planeB.x*factorB, planeB.y*factorB, planeB.z*factorB, planeB.w*factorB);
		j = tempPlaneA.y - tempPlaneB.y;
		k = tempPlaneA.z - tempPlaneB.z;
		l = tempPlaneA.w - tempPlaneB.w;
		//System.out.println(j+", "+k+", "+l);
		z = (l-j)/k;
		Vector3f planeC = new Vector3f(planeA.x,planeA.y+planeA.z*z,planeA.w);
		Vector3f planeD = new Vector3f(planeB.x,planeB.y+planeB.z*z,planeB.w);
		factor = planeC.x*planeD.x;
		float factorC = factor/planeC.x;
		float factorD = factor/planeD.x;
		if (factor == 0){
			factorC = 1;
			factorD = 1;
		}
		tempPlaneC = new Vector3f(planeC.x*factorC, planeC.y*factorC, planeC.z*factorC);
		tempPlaneD = new Vector3f(planeD.x*factorD, planeD.y*factorD, planeD.z*factorD);
		tempPlaneX = new Vector3f(tempPlaneC.x-tempPlaneD.x, tempPlaneC.y-tempPlaneD.y, tempPlaneC.z-tempPlaneD.z);
		y = tempPlaneX.z/tempPlaneX.y;
		x = (planeC.z - (planeC.y*y))/planeC.x;
		z = (planeA.w -(planeA.x*x) - (planeA.y*y))/planeA.z;
		return new Vector3f(x,y,z);
	}

	private Vector3f comparePlanes(Vector3f plane10, Vector3f plane11, Vector3f plane12, Vector3f plane20 , Vector3f plane21 , Vector3f plane22){
		//System.out.println("10."+plane10);
		//System.out.println("11."+plane11);
		//System.out.println("12."+plane12);
		///System.out.println("20."+plane20);
		//System.out.println("21."+plane21);
		//System.out.println("22."+plane22);
		Vector3f a1 = new Vector3f((plane11.x-plane10.x),(plane11.y-plane10.y),(plane11.z-plane10.z));
		Vector3f b1 = new Vector3f((plane12.x-plane10.x),(plane12.y-plane10.y),(plane12.z-plane10.z));
		Vector3f crossPorduct = crossProduct(a1, b1);
		Vector4f planeA = formatPlane(plane10, crossPorduct);
		Vector3f a2 = new Vector3f((plane21.x-plane20.x),(plane21.y-plane20.y),(plane21.z-plane20.z));
		Vector3f b2 = new Vector3f((plane22.x-plane20.x),(plane22.y-plane20.y),(plane22.z-plane20.z));
		//System.out.println(crossPorduct);
		crossPorduct = crossProduct(a2, b2);
		Vector4f planeB = formatPlane(plane20, crossPorduct);
		//System.out.println(crossPorduct);
		return intersection(planeA,planeB);
	}


}
