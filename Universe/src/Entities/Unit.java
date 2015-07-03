package Entities;

import org.lwjgl.util.vector.Vector3f;
import org.lwjgl.util.vector.Vector4f;

import renderEngine.TexturedModel;

public class Unit {
	private int unitID;
	private TexturedModel model;
	private Vector3f position;
	private Vector3f angle;
	private Vector3f angularSpeed;
	private Vector3f directionalSpeed;
	private float scale;
	private float strength;
	private float energy = 100000;
	private float electronegativity;
	private float mass = 100;
	private float connectionStrenth;
	private float connectionFlexibility;
	private Vector4f colour = new Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
	private float dampner;
	private float reflectivity;
	private boolean isLight;
	private boolean useFakeLighting;
	
	public Unit(TexturedModel model, Vector3f position, Vector3f angle, Vector3f directionalSpeed, Vector3f angularSpeed, float scale){
		this.model = model;
		this.position = position;
		this.angle = angle;
		this.directionalSpeed = directionalSpeed;
		this.angularSpeed = angularSpeed;
		this.scale = scale;
	}
	public Vector3f rotateUnit(Vector3f point){
		Vector3f result = new Vector3f(0,0,0);
		//TRANSLATE TO ORIGIN
		float x1 = (point.x - position.x);
		float y1 = (point.y - position.y);

		//APPLY ROTATION
		float temp_x1 = (float) (x1 * Math.cos(Math.toRadians(angle.z)) - y1 * Math.sin(Math.toRadians(angle.z)));
		float temp_y1 = (float) (x1 * Math.sin(Math.toRadians(angle.z)) + y1 * Math.cos(Math.toRadians(angle.z)));

		//TRANSLATE BACK
		point.x = (float) (temp_x1 + position.x);
		point.y = (float) (temp_y1 + position.y);

		//TRANSLATE TO ORIGIN
		x1 = (point.x - position.x);
		float z1 = (point.z - position.z);

		//APPLY ROTATION
		temp_x1 = (float) (x1 * Math.cos(Math.toRadians(angle.y)) - z1 * Math.sin(Math.toRadians(angle.y)));
		float temp_z1 = (float) (x1 * Math.sin(Math.toRadians(angle.y)) + z1 * Math.cos(Math.toRadians(angle.y)));

		//TRANSLATE BACK
		point.x = temp_x1 + position.x;
		point.z = temp_z1 + position.z;
		
		//TRANSLATE TO ORIGIN
		y1 = (point.y - position.y);
		z1 = (point.z - position.z);

		//APPLY ROTATION
		temp_y1 = (float) (y1 * Math.cos(Math.toRadians(angle.x)) - z1 * Math.sin(Math.toRadians(angle.x)));
		temp_z1 = (float) (y1 * Math.sin(Math.toRadians(angle.x)) + z1 * Math.cos(Math.toRadians(angle.x)));

		//TRANSLATE BACK
		point.y = temp_y1 + position.y;
		point.z = temp_z1 + position.z;
		
		result = new Vector3f(point.x,point.y,point.z);
		return result;
	}
	
	public Vector3f[] getUnitPoints(){
		Vector3f FrontPoint = new Vector3f(
				(float) (scale * Math.cos(Math.toRadians(90))) + position.x, 
				(float) (scale * Math.cos(Math.toRadians(135)))+ position.y, 
				(float) (scale * Math.cos(Math.toRadians(180)))+ position.z);
		FrontPoint = rotateUnit(FrontPoint);
		Vector3f RightPoint = new Vector3f(
				(float) (scale * Math.cos(Math.toRadians(30)))+ position.x, 
				(float) (scale * Math.cos(Math.toRadians(135)))+ position.y, 
				(float) (scale * Math.cos(Math.toRadians(60)))+ position.z);
		RightPoint = rotateUnit(RightPoint);
		Vector3f TopPoint = new Vector3f(
				(float) (scale * Math.cos(Math.toRadians(90)))+ position.x, 
				(float) (scale * Math.cos(Math.toRadians(45)))+ position.y, 
				(float) (scale * Math.cos(Math.toRadians(90)))+ position.z);
		TopPoint = rotateUnit(TopPoint);
		Vector3f LeftPoint = new Vector3f(
				(float) (scale * Math.cos(Math.toRadians(150)))+ position.x, 
				(float) (scale * Math.cos(Math.toRadians(135)))+ position.y, 
				(float) (scale * Math.cos(Math.toRadians(60)))+ position.z);
		LeftPoint = rotateUnit(LeftPoint);
		Vector3f[] results = new Vector3f[4];
		results[0] = FrontPoint;
		results[1] = RightPoint;
		results[2] = TopPoint;
		results[3] = LeftPoint;
		return results;
	}
	
	public void update(){
		increaseRotation(angularSpeed.x, angularSpeed.y, angularSpeed.z);
		increasePosition(directionalSpeed.x, directionalSpeed.y, directionalSpeed.z);
	}
	
	public void increasePosition(float dx, float dy, float dz){
		this.position.x += dx;
		this.position.y += dy;
		this.position.z += dz;
	}
	
	public void increaseRotation(float dx, float dy, float dz){
		this.angle.x += dx;
		this.angle.y += dy;
		this.angle.z += dz;
	}
	
	public void increaseDirectionalSpeed(float dx, float dy, float dz){
		this.directionalSpeed.x += dx;
		this.directionalSpeed.y += dy;
		this.directionalSpeed.z += dz;
	}
	
	public void increaseAngularSpeed(float dx, float dy, float dz){
		this.angularSpeed.x += dx;
		this.angularSpeed.y += dy;
		this.angularSpeed.z += dz;
	}
	
	
	
	public boolean isLight() {
		return isLight;
	}

	public void setLight(boolean isLight) {
		this.isLight = isLight;
	}

	public boolean isUseFakeLighting() {
		return useFakeLighting;
	}

	public void setUseFakeLighting(boolean useFakeLighting) {
		this.useFakeLighting = useFakeLighting;
	}

	public float getScale() {
		return scale;
	}

	public void setScale(float scale) {
		this.scale = scale;
	}
	
	public int getUnitID() {
		return unitID;
	}

	public void setUnitID(int unitID) {
		this.unitID = unitID;
	}

	public float getStrength() {
		return strength;
	}

	public void setStrength(float strength) {
		this.strength = strength;
	}

	public Vector3f getAngularSpeed() {
		return angularSpeed;
	}

	public void setAngularSpeed(Vector3f angularSpeed) {
		this.angularSpeed = angularSpeed;
	}

	public Vector3f getDirectionalSpeed() {
		return directionalSpeed;
	}

	public void setDirectionalSpeed(Vector3f directionalSpeed) {
		this.directionalSpeed = directionalSpeed;
	}

	public TexturedModel getModel() {
		return model;
	}

	public void setModel(TexturedModel model) {
		this.model = model;
	}

	public Vector3f getPosition() {
		return position;
	}

	public void setPosition(Vector3f position) {
		this.position = position;
	}

	public Vector3f getAngle() {
		return angle;
	}

	public void setAngle(Vector3f angle) {
		this.angle = angle;
	}


	public float getEnergy() {
		return energy;
	}

	public void setEnergy(float energy) {
		this.energy = energy;
	}

	public float getElectronegativity() {
		return electronegativity;
	}

	public void setElectronegativity(float electronegativity) {
		this.electronegativity = electronegativity;
	}

	public float getMass() {
		return mass;
	}

	public void setMass(float mass) {
		this.mass = mass;
	}

	public float getConnectionStrenth() {
		return connectionStrenth;
	}

	public void setConnectionStrenth(float connectionStrenth) {
		this.connectionStrenth = connectionStrenth;
	}

	public float getConnectionFlexibility() {
		return connectionFlexibility;
	}

	public void setConnectionFlexibility(float connectionFlexibility) {
		this.connectionFlexibility = connectionFlexibility;
	}

	public float getDampner() {
		return dampner;
	}

	public void setDampner(float dampner) {
		this.dampner = dampner;
	}

	public float getReflectivity() {
		return reflectivity;
	}

	public void setReflectivity(float reflectivity) {
		this.reflectivity = reflectivity;
	}

	/**
	 * @return the colour
	 */
	public Vector4f getColour() {
		return colour;
	}

	/**
	 * @param colour the colour to set
	 */
	public void setColour(Vector4f colour) {
		this.colour = colour;
	}
	
	

}
