package less.gui.model;

import javafx.beans.property.SimpleStringProperty;

public class PositionXY {
	
	private SimpleStringProperty pos_x;
	private SimpleStringProperty pos_y;
	private SimpleStringProperty pos_z;// the z position of object
	
	private SimpleStringProperty extra_props; //extra field, current, it contains rotation
	
	public PositionXY(String x, String y){
		this.pos_x = new SimpleStringProperty(x);
		this.pos_y = new SimpleStringProperty(y);
		this.pos_z = new SimpleStringProperty("0");// rotation
		this.extra_props = new SimpleStringProperty("0 ");// by default, it adds a rotation degree
	}
	
	public PositionXY(String x, String y, String z){
		this.pos_x = new SimpleStringProperty(x);
		this.pos_y = new SimpleStringProperty(y);
		this.pos_z = new SimpleStringProperty(z);
		this.extra_props = new SimpleStringProperty("0 ");
	}
	
	public PositionXY(String x, String y, String z, String props){
		this.pos_x = new SimpleStringProperty(x);
		this.pos_y = new SimpleStringProperty(y);
		this.pos_z = new SimpleStringProperty(z);
		this.extra_props = new SimpleStringProperty(props);
	}
	
	public String getPos_x() {
		return pos_x.get();
	}
	public void setPos_x(String pos_x) {
		this.pos_x.set(pos_x);
	}
	public String getPos_y() {
		return pos_y.get();
	}
	public void setPos_y(String pos_y) {
		this.pos_y.set(pos_y);
	}
	
	public String getPos_z() {
		return pos_z.get();
	}
	public void setPos_z(String pos_z) {
		this.pos_z.set(pos_z);
	}
	
	public String getExtra_props(){
		return extra_props.get();
	}
	
	public void setExtra_props(String extra_props){
		this.extra_props.set(extra_props);
	}
	
	
}
