package less.gui.utils;

import java.util.ArrayList;

import javafx.scene.paint.Color;

/**
 * 
 * @author Jim
 *
 */
public class Const {
	public static String LESS_TITLE = "LESS-";
	public static String LESS_LAUNCH_SCRIPT= "less.py";
	public static String LESS_LAUNCH_SCRIPT_PYC= "less.pyc";
	public static String LESS_LAUNCH_PATH = "E:\\Coding\\LessProject\\python\\lesspy\\";
	public static String LESS_CONST_JSON_NAME="const.conf";
	public static String LESS_BATCH_JSON_NAME="batch.conf";
	public static String LESS_DOT_LESS_FILE = ".less";
	public static String LESS_CONFIG_PROPERTIES_FILE="cfg.properties";
	
	public static String LESS_DEFAULT_INPUT_FILE = "default.conf";
	
	public static String LESS_TERRAIN_MESH = "MESH";
	public static String LESS_TERRAIN_PLANE = "PLANE";
	public static String LESS_TERRAIN_RASTER = "RASTER";
	public static String LESS_TERRAIN_LANDCOVER_FILE = "landcover.txt";
	
	public static String LESS_SENSOR_TYPE_ORTH = "orthographic";
	public static String LESS_SENSOR_TYPE_PER = "perspective";
	public static String LESS_SENSOR_TYPE_PT = "PhotonTracing";
	
	public static String LESS_OBJECTS_FILE_NAME = "objects.txt"; //for store objects of tree models.
	public static String LESS_INSTANCE_FILE_NAME = "instances.txt"; //for store the position of tree models.
	public static String LESS_HIDE_OBJECT_FILE_NAME = "hide_object.txt"; //for hide some objects
	public static String LESS_OBJECTS_ALTITUDE_3D = "object_pos_3dView.txt"; //
	public static String LESS_OBJECTS_BOUNDINGBOX_FILE = "objects_boundingbox.txt";

	public static int LESS_TREE_POS_DRAW_MAX_NUM = 1000000;// 
	
	public static String LAST_OPENED_SIM = "simu";
	public static String LAST_OPNED_CHOOSE_OBJ = "chooseObj";
	public static String LAST_OPNED_CHOOSE_LANDCOVER = "chooselc";
	public static String LAST_OPNED_CHOOSE_BACKIMG = "choosebackimg";
	public static String LAST_OPNED_CHOOSE_POLYGON = "choosePolygon";
	public static String LAST_OPNED_DEF_FILE_PATH = "chooseDef";
	public static String LAST_OPNED_BATCH_FILE_PATH = "chooseBatch";
	public static String LAST_OPNED_CHM_FILE_PATH = "chooseCHM";
	
	public static String LESS_RT_NAME_WINDOWS = "lessrt.exe";
	public static String LESS_RT_NAME_LINUX = "lessrt";
	
	
	// LESS Mode
	public static String LESS_MODE= "development";
	public static String LESS_VERSION = "V1.6";
	
	public static boolean LESS_OUT_ALL = false;
	
	public static boolean LESS_HIDE_NOT_IMPLEMENTED = true;

	public static String LESS_TREE_POS_DIS_UNIFORM = "Uniform";
	public static String LESS_TREE_POS_DIS_POISSON = "Poisson";
	
	
	//database
	public static String LESS_DBNAME_LambertianDB = "LambertianDB.db";
	
	//image format
	public static String LESS_OUT_IMAGEFORMAT_SPECTRUM = "spectrum";
	public static String LESS_OUT_IMAGEFORMAT_RGB = "rgb";	
	
	//mask
	public static String LESS_SCRIPTS_MASK = "java_domask";	
	
	//Utility
	public static String LESS_UTILITY_SCRIPT_DEF2OBJ = "def2obj";
	public static String LESS_UTILITY_SCRIPT_LCTYPES = "extractLandcoverType";
	public static String LESS_UTILITY_SCRIPT_RASTER2OBJ = "Raster2obj";
	
	public static String LESS_PYTHON_STARTUP_SCRIPT = "startup_scripts";
	
	public static String LESS_BATCH_TOOL_TITLE="Batch Tool";
	public static String LESS_BATCH_SCRIPT_NAME="BatchGenAndRun";
	public static String LESS_SCRIPT_POST_PROCESSING="PostProcessing";
	public static String LESS_SCRIPT_TREE_DETECTION="Singletree_detection_from_CHM";
	
	//optical
	public static String LESS_DEFAULT_OPTICAL1 = "birch_branch";
	public static String LESS_DEFAULT_OPTICAL2 = "dark_soil_mollisol";
	public static String LESS_DEFAULT_OPTICAL3 = "birch_leaf_green";
	public static String LESS_DEFAULT_TEMPERATURE = "T300";
	public static Color LESS_DEFAULT_BRANCH_COLOR = Color.web("#b34d1a");
	
	
	//cluster
	public static String LESS_SERVER_EXE = "lesssrv.exe";
	public static String LESS_SERVER_TXT_FILE = "server.txt";
	
	public static String LESS_HELP_ONLINE_URL = "http://ramm.bnu.edu.cn/projects/less/Documentation/";
	
	//LAI
	public static String LESS_LAI_OUTPUT_FILE = "LAI.txt";
	
}
