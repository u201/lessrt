package less.gui.utils;

import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.TableView;
import less.gui.model.FacetOptical;

public class DBReader {
	
	private Connection c = null;
    private Statement stmt = null;
	
	public String getLambertDBPath(){
		String dbpath;
		if(Const.LESS_MODE.equals("development")){
			dbpath= Paths.get(Const.LESS_LAUNCH_PATH,"SpectralDB",Const.LESS_DBNAME_LambertianDB).toString();
		}else{
			dbpath = Paths.get(System.getProperty("user.dir"),"bin","scripts","Lesspy",
					"SpectralDB",Const.LESS_DBNAME_LambertianDB).toString();
		}
		return dbpath;
	}
	
	public void connect(){	    
	    try {
	    	Class.forName("org.sqlite.JDBC");
			c = DriverManager.getConnection("jdbc:sqlite:"+getLambertDBPath());
			stmt = c.createStatement();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public ObservableList<String> getTableList(){
		ObservableList<String> tableList = FXCollections.observableArrayList();
		try {
			ResultSet rs = stmt.executeQuery( "SELECT name FROM sqlite_sequence;");
			while ( rs.next() ) {
				tableList.add(rs.getString("name"));
			}
			rs.close();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return tableList;
	}
	
	/**
	 * 根据光谱名字，选择波段的反射率或者透过率
	 * @param opticalName 光谱名字
	 * @return
	 */
	public ArrayList<String> getOpticalByName(String opticalName, String wavelength_and_bandwidth){
		ArrayList<String> opticalArr = new ArrayList<String>();
		String[] wAbarr = wavelength_and_bandwidth.split(",");
		String rfString="",bfString = "", transString="";
		for(int i=0;i<wAbarr.length;i++){
			String[] wb = wAbarr[i].split(":");
			double center_w = Double.parseDouble(wb[0]);
			double bandwidth = Double.parseDouble(wb[1]);
//			if(bandwidth<1)
//				bandwidth = 1;
			double left = center_w - 0.5*bandwidth;
			double right = center_w + 0.5*bandwidth;
			double fr=0,br=0,trans=0,num=0;
			try {
				String sql = "SELECT wavelength, front_ref,back_ref,transmittance FROM "+opticalName+" where wavelength>="+
						left + " and wavelength<="+ right + ";";
				ResultSet rs = stmt.executeQuery(sql);
				while ( rs.next() ) {
					num++;
					double front_ref = rs.getDouble("front_ref");
					fr += front_ref;
					double back_ref = rs.getDouble("back_ref");
					br += back_ref;
					double transmittance = rs.getDouble("transmittance");
					trans += transmittance;
				}
				if(num!=0){
					fr /= num;
					br /= num;
					trans /= num;
				}else{
					//当没有找到合适的波段是，则进行插值.
					sql = "SELECT wavelength, front_ref,back_ref,transmittance FROM "+opticalName+" order by abs(wavelength-"+center_w
							+ ") LIMIT 2";
					rs = stmt.executeQuery(sql);
					rs.next();
					double front_ref1 = rs.getDouble("front_ref");
					double back_ref1 = rs.getDouble("back_ref");
					double transmittance1 = rs.getDouble("transmittance");
					rs.next();
					double front_ref2 = rs.getDouble("front_ref");
					double back_ref2 = rs.getDouble("back_ref");
					double transmittance2 = rs.getDouble("transmittance");
					fr = 0.5*(front_ref1+front_ref2); //应该用线性插值，不过现在就简单平均吧
					br = 0.5*(back_ref1+back_ref2);
					trans = 0.5 * (transmittance1+transmittance2);
				}
				
				rfString += String.format("%.5f", fr) +",";
				bfString += String.format("%.5f", br)+",";
				transString += String.format("%.5f", trans)+",";
				rs.close();
				
			} catch (SQLException e) {
				e.printStackTrace();
			}
		}
		opticalArr.add(rfString.substring(0, rfString.length()-1));
		opticalArr.add(bfString.substring(0, bfString.length()-1));
		opticalArr.add(transString.substring(0, transString.length()-1));
		return opticalArr;
	}
	
	public FacetOptical getOpticalTableElementByName(String opticalName, String wavelength_and_bandwidth){
		ArrayList<String> optical = this.getOpticalByName(opticalName, wavelength_and_bandwidth);
		return new FacetOptical(opticalName, optical.get(0), optical.get(1), optical.get(2));
	}
	
	/**
	 * 根据光谱库中的定义的进行更新，如果改变波段的情况下
	 * @param opticalData
	 */
	public void refreshOpticalDB(ObservableList<FacetOptical> opticalData,String waveLength_and_bandwidth){
		ObservableList<String> tables = this.getTableList();
		ObservableList<FacetOptical> opticalDatatmp = FXCollections.observableArrayList();
		for(int i=0;i<opticalData.size();i++){
			String opticalName = opticalData.get(i).getOpticalName();
			if(tables.contains(opticalName)){
				ArrayList<String> opticalArr = this.getOpticalByName(opticalName, waveLength_and_bandwidth);
				opticalDatatmp.add(new FacetOptical(
						opticalName,
						opticalArr.get(0),
						opticalArr.get(1),
						opticalArr.get(2)
						));
			}
			else{
				opticalDatatmp.add(opticalData.get(i));
			}
		}
		opticalData.clear();
		opticalData.addAll(opticalDatatmp);
	}
}
