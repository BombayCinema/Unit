package Entities;

import org.lwjgl.input.Keyboard;
import org.lwjgl.input.Mouse;
import org.lwjgl.util.vector.Vector3f;

import renderEngine.DisplayManager;
import renderEngine.TexturedModel;

public class Player extends Unit {

	private static final float X_SPEED = 1;
	private static final float Y_SPEED = 1;
	private static final float Z_SPEED = 1;

	private float currentXSpeed = 0;
	private float currentYSpeed = 0;
	private float currentZSpeed = 0;
	private float currentRotation = 0;


	public Player(TexturedModel model, Vector3f position, float rotX,
			float rotY, float rotZ, float scale) {
		super(model, position, new Vector3f(rotX, rotY, rotZ),  new Vector3f(0, 0, 0), new Vector3f(0, 0, 0), scale);
	}

	public void move(){
		checkInputs();
		float speed;
		float tempAngle = (float) (Math.atan2(currentXSpeed, currentZSpeed) * 180 / Math.PI);
		if (currentXSpeed !=0 || currentZSpeed !=0){
			speed = 1*DisplayManager.getFrameTimeSeconds();
		} else {
			speed = 0;
		}
		float dx = (float) (speed * Math.sin(Math.toRadians(getAngle().y+tempAngle)));
		float dz = (float) (speed * Math.cos(Math.toRadians(getAngle().y+tempAngle)));
		increaseDirectionalSpeed(dx, currentYSpeed*DisplayManager.getFrameTimeSeconds(), dz);
		increaseRotation(0, currentRotation*DisplayManager.getFrameTimeSeconds(), 0);
	}



	private void checkInputs(){
		if (Keyboard.isKeyDown(Keyboard.KEY_W)){
			this.currentZSpeed = -Z_SPEED;
		} else if (Keyboard.isKeyDown(Keyboard.KEY_S)){
			this.currentZSpeed = Z_SPEED;
		} else {
			this.currentZSpeed = 0;
		}

		if (Keyboard.isKeyDown(Keyboard.KEY_D)){
			this.currentXSpeed = X_SPEED;
		} else if (Keyboard.isKeyDown(Keyboard.KEY_A)){
			this.currentXSpeed = -X_SPEED;
		} else {
			this.currentXSpeed = 0;
		}

		if (Keyboard.isKeyDown(Keyboard.KEY_LSHIFT)){
			this.currentYSpeed = Y_SPEED;
		} else if (Keyboard.isKeyDown(Keyboard.KEY_LCONTROL)){
			this.currentYSpeed = -Y_SPEED;
		} else {
			this.currentYSpeed = 0;
		}

		if(Mouse.isButtonDown(0)) {
			float pitchChange = Mouse.getDX() * 1f;
			currentRotation -= pitchChange;
		} else {
			currentRotation = 0;
		}
	}
}
