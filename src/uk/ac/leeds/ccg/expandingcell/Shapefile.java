package uk.ac.leeds.ccg.expandingcell;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFactorySpi;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FeatureSource;
import org.geotools.data.FeatureStore;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureCollections;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

/**
 * Represents a shape file (used by ArcGIS). Basically a convenience function to simplify the
 * tasks of reading/writing shapefiles (using GeoTools). Most of the shapefile reading/writing/editing 
 * code has been copied or adapted from the GeoTools tutorials which come with the GeoToolssoftware.
 * <p>
 * Typical usage (reads a shapefile, get an attribute, adds a new attribute column, then writes it):<br>
 * <code>
 * Shapefile s = new Shapefile("ObjectIDs", "myShapefile.shp");<br>
 * s.readShapefile();<br>
 * int anAttributeValue = s.getAttribute(featureID, "anAttribute", Integer.class);<br>
 * s.addAttribute("newAttribute", Double.class, valuesTable)
 * s.writeShapefile("newShapefile");
 * </code>
 * 
 * @author Nick Malleson
 */
public class Shapefile implements Cloneable {

	/** The geometries of the objects in this Shapefile, stored along with object ids.*/
	private Hashtable<Integer, Geometry> geometries;

	/** The features associated with each object in this Shapefile, stored along with object ids. */
	private Hashtable<Integer, SimpleFeature> features;

	/** The actual file which this shapefile has been created from */
	private File file;

	/** The FeatureSource for this shapefile, required as a template for writing shapefiles */
	private FeatureSource<SimpleFeatureType, SimpleFeature> featureSource;

	/** Optional column used to store object ids (by default will use FID). */
	private String idColumn = null;

	/** When creating new shapefiles this will sometimes be set */
	private DefaultGeographicCRS crs;

	/**
	 * Create a new Shapefile object.
	 * @param idColumn
	 * @param file
	 */
	public Shapefile(String idColumn, File file) {
		this.geometries = new Hashtable<Integer, Geometry>();
		this.features = new Hashtable<Integer, SimpleFeature>();
		this.idColumn = idColumn;
		this.file = file;
	}

	/**
	 * Create a new Shapefile object, using the default FID as the ID for each object.
	 * @param file
	 */
	public Shapefile(File file) {
		this.geometries = new Hashtable<Integer, Geometry>();
		this.features = new Hashtable<Integer, SimpleFeature>();
		this.file = file;
	}


	/** Default constructor is used by the <code>copy</code> function. */
	private Shapefile() { 
	}

	/**
	 * Get the value of an attribute.
	 * @param featureID The ID of the object to get the attribute value from.
	 * @param attributeName The name of the attribute.
	 * @param clazz The class of the obect to be returned.
	 * @return The attribute as an Object or 
	 */
	@SuppressWarnings("unchecked")
	public <T> T getAttribute(int featureID, String attributeName, Class<T> clazz) {
		if (this.checkFeatures()) { // Check there are some features in this Shapefile
			SimpleFeature feature = this.features.get(featureID);
			if (feature == null) {
				System.err.println("Shapefile.getAttribute() error: no feature with id: "+featureID+
				" in this Shapefile");
				return null;
			}
			Object attribute = feature.getAttribute(attributeName);
			if (attribute == null) {
				System.err.println("Shapefile.getAttribute() error: no colum named: "+attributeName);
				return null;
			}
			T a = null; // The attribute cast to it's correct class.
			if (attribute.getClass().isAssignableFrom(clazz)) {
				a = (T) attribute;
			}
			else {
				System.err.println("The attribute "+attributeName+" cannot be assigned to a "+
						clazz.getName()+", it is a: "+attribute.getClass().getName());
			}
			return a;
		}
		else { // No features in this Shapefile
			System.err.println("Shapefile.getAttribute() error: no features in this shapefile");
			return null;
		}
	}

	/**
	 * Add a new attribute to this Shapefile.
	 * @param attributeName The name of the attribute to add.
	 * @param valuesClass The class of the values to be added.
	 * @param attributeValues A map of the IDs of each feature and the associated value of the
	 * new attribute
	 * @param T The type of attribute to add.
	 * @return True if the operation was successful, false otherwise.
	 */
	public <T> boolean addAttribute(String attributeName, Class<T> valuesClass, Map<Integer, T> attributeValues) {
		if (!this.checkFeatures()) {
			System.err.println("Shapefile.addAttribute() error: no features in this Shapefile ("+this.file.getAbsolutePath()+")");
			return false;
		}
		// Do some checks on the input
		boolean error = false;
		// Check nothing is null
		if (attributeName == null ){
			System.err.println("Shapefile.addAttribute() error: the input attribute name is null.");
			error = true;
		}
		if (valuesClass == null ){
			System.err.println("Shapefile.addAttribute() error: the input values class is null.");
			error = true;
		}
		if (attributeValues == null ) {
			System.err.println("Shapefile.addAttribute() error: the input attribute values are null");
			error = true;
		}
		// Check the input values are ok (ID's match exactly)
		if (attributeValues.size() != this.features.size()) {
			System.err.println("Shapefile.addAttribute() error: there are "+attributeValues.size()+" input values "+
					"and "+this.features.size()+" features in this shapefile. Cannot add attributes unless these " +
			"numbers are the same.");
			error = true;
		}
		List<Integer> badIDs = new ArrayList<Integer>();
		for (Integer i:attributeValues.keySet())
			if (!this.features.containsKey(i))
				badIDs.add(i);
		if (badIDs.size()>0) {
			System.err.println("Shapefile.addAttribute() error: could not find features associated " +
					"with the following IDs: "+badIDs.toString());
			error = true;
		} // if badIDs
		if (error) { // If there has been an error return false
			System.err.println("Shapefile.addAttribute(): there are problems with the inputs to this function " +
			"(see previous error messages). Cannot continue adding this new attribute.");
			return false;
		}

		/* Method works be re-building all features with the new attributed added */
		// Create a feature type builder to create new features with the added attribute
		SimpleFeatureTypeBuilder featureTypeBuilder = new SimpleFeatureTypeBuilder();
		// Use first feature to initialise the feature builder
		featureTypeBuilder.init(features.values().iterator().next().getFeatureType());
		// Add the new attribute
		featureTypeBuilder.add(attributeName, valuesClass);
		// Create a feature builder to create the new features
		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureTypeBuilder.buildFeatureType());
		// Iterate over all existing features, creating new ones with the added attribute
		Map<Integer, SimpleFeature> newFeatures = new Hashtable<Integer, SimpleFeature>();
		for (Integer id:this.features.keySet()) {
			SimpleFeature newFeature = featureBuilder.buildFeature(String.valueOf(id));
			SimpleFeature existingFeature = this.features.get(id);
			// Add all existing attributes to the new feature
			for (int i = 0; i < existingFeature.getAttributeCount(); i++) {
				newFeature.setAttribute(i, existingFeature.getAttribute(i));
			}
			// Add the new attribute to the new feature
			newFeature.setAttribute(attributeName, attributeValues.get(id));
			// Replace the existing feature with the new one
			newFeatures.put(id, newFeature);
		}
		// Finally replace all old features with the old ones
		for (Integer id:this.features.keySet()) {
			this.features.put(id, newFeatures.get(id));
		}
		return true;
	}

	/** Make sure there are some features in this Shapefile, return fals if not. */
	private boolean checkFeatures() {
		if (this.features==null || this.features.size()==0) 
			return false;
		else
			return true;
	}

	/**
	 * Read in all objects stored in the shapefile, adding them to this Shapefile object.
	 * @return true if everything was read in successfully, false otherwise.
	 */
	public boolean readShapefile() {

		this.clear(); // Clear all objects from this Shapefile

		// Connection to the shapefile
		Map<String, Serializable> connectParameters = new HashMap<String, Serializable>();

		try {
			connectParameters.put("url", this.file.toURI().toURL());
			connectParameters.put("create spatial index", true);
			DataStore dataStore = DataStoreFinder.getDataStore(connectParameters);

			// we are now connected
			String[] typeNames = dataStore.getTypeNames();
			String typeName = typeNames[0];

			this.featureSource = dataStore.getFeatureSource(typeName); 
			FeatureCollection<SimpleFeatureType, SimpleFeature> collection = featureSource.getFeatures();
			FeatureIterator<SimpleFeature> iterator = collection.features();

			try {
				while (iterator.hasNext()) {
					SimpleFeature feature = iterator.next();
					// Set the feature's ID
					int id = 0;
					try {
						if (this.idColumn!=null) // User has supplied an ID column
							id = (Integer)feature.getAttribute(this.idColumn);
						else // Use the FID colum to uniquely identify each object (stored as 'Name.FID' so split on '.')
							id = Integer.parseInt(feature.getID().split("\\.")[1]);
					}
					catch (ClassCastException e) {
						System.err.println("Shapfile.readObjects() error: cannot read integer ids from the " +
								"column: "+this.idColumn+". Check this column stores unique integer IDs.");
						return false;
					}
					catch (NullPointerException e) {
						System.err.println("Shapfile.readObjects() error: cannot read integer ids from the" +
								"column: "+this.idColumn+". Check this column stores unique integer IDs.");
						return false;
					}
					catch (NumberFormatException e) {
						System.err.println("Shapfile.readObjects() error: cannot cast this feature's " +
								"ID to an integer: "+feature.getID().toString());
						return false;
					}

					if (features.containsKey(id)) {
						System.err.println("Shapefile.readObjects() error: this feature's ID ("+id
								+")is not unique ");
						return false;
					}
					features.put(id, feature);
					Geometry geometry = (Geometry) feature.getDefaultGeometry();
					geometries.put(id, geometry);                    
				} // while iterator.hasNext()
			} finally {
				if (iterator != null) {
					iterator.close();
				}
			}
		} catch (MalformedURLException e) {
			e.printStackTrace();
			return false;
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		return true;
	} // readShapefile

	/**
	 * Creates point features for this shapefile. This code has been adapted from the GeoTools 2.5.1 tutorial
	 * documentation ("CSV2SHP Lab" chapter).
	 * @param xcoords A list of x coordinates.
	 * @param ycoords A list of y coordinates.
	 * @param ids A list of id's for each feature
	 * @param idColumnName The name of the column which holds the ids
	 * @param crs The coordinate reference system for the new features (if null, WGS84 will be used).
	 * @return True if the operation was successful, false otherwise.
	 */
	public boolean createPointFeatures(List<Double> xcoords, List<Double> ycoords, List<Integer> ids, 
			String idColumnName, DefaultGeographicCRS crs) {

		// Check this shapefile doesn't already have features in it
		if (this.features.size() > 0 ) {
			System.err.println("Shapefile.createPointFeatures() error: this Shapefile already has features in it, " +
			"it must be cleared before new features can be added");
			return false;
		}
		// Check all lists are same length
		if (xcoords.size() != ycoords.size() || xcoords.size() != ids.size()) {
			System.err.println("Shapefile.createPointFeatures() error, lengths of lists do not match, they must "+
			"all be the same size");
			return false;
		}	
		// Remember the CRS for later (if writeShapefile() is called).
		this.crs = crs;
		//		try {
		// Create the features and geometries
		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName( "Location" );
		builder.setCRS( this.crs == null ? DefaultGeographicCRS.WGS84 : this.crs );
		builder.add( "Location", Point.class );
		builder.add( "ID", Integer.class);
		SimpleFeatureType TYPE = builder.buildFeatureType();
		// (Old way of creating feature type, not as good because cannot set CRS)
		//			SimpleFeatureType TYPE = DataUtilities.createType("Location", "location:Point,"+idColumnName+":Integer");
		GeometryFactory geomFac = JTSFactoryFinder.getGeometryFactory(null);
		for (int i=0; i<xcoords.size(); i++) {
			int id = ids.get(i);
			Point point = geomFac.createPoint( new Coordinate(xcoords.get(i),ycoords.get(i)));
			SimpleFeature feature = SimpleFeatureBuilder.build( TYPE, new Object[]{point, id}, null );
			//				collection.add( feature );
			this.features.put(id, feature);
			this.geometries.put(id, point);
		}

		//		} catch (SchemaException e) {
		//			System.err.println("Shapefile.createPointFeatures error: "+(e.getMessage()!=null ? e.getMessage() : ""));
		//			e.printStackTrace();
		//			return false;
		//		}
		return true;
	}

	/**
	 * Create features for the shapefile using the given geometries.
	 * @param geometries The geometries to create the features from
	 * @param ids The id's of each feature
	 * @param idColumnName The column name to hold the ids
	 * @param crs The coordinate reference system to use (if null, WGS1984 will be used).
	 * @return True if successful, false otherwise.
	 */
	public boolean createFeatures(List<Geometry> geometries, List<Integer> ids, String idColumnName,
			DefaultGeographicCRS crs ) {
		// Check this shapefile doesn't already have features in it
		if (this.features.size() > 0 ) {
			System.err.println("Shapefile.createFeatures() error: this Shapefile already has features in it, " +
			"it must be cleared before new features can be added");
			return false;
		}
		// Check the lists aren't empty
		if (geometries == null || geometries.size() == 0) {
			System.err.println("Shapefile.createFeatures() error, the list of geometries is null or empty.");
			return false;
		}
		// Check lists are same length
		if (geometries.size() != ids.size()) {
			System.err.println("Shapefile.createFeatures() error, the length of the geometries list is " +
			"a differnet length to the list of ids");
			return false;
		}
		// Check the geometries are all the same
		Geometry testGeom = null;
		for (int i=0; i<geometries.size(); i++) {
			if (i==0) {
				testGeom = geometries.get(0);
			}
			else {
				if (!testGeom.getClass().equals(geometries.get(i).getClass())) {
					System.err.println("Shapefile.createFeatures() error, not all geometries are the same, have " +
							"found a "+testGeom.getClass().getName()+" and a "+geometries.get(i).getClass().getName());
					return false;
				}
			}
		} // for geometries
		// Remember the CRS for later (if writeShapefile() is called).
		this.crs = crs;
		// Create the features and geometries
		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName( "Location" );
		builder.setCRS( this.crs == null ? DefaultGeographicCRS.WGS84 : this.crs );
		builder.add( "Location", testGeom.getClass());
		builder.add( "ID", Integer.class);
		SimpleFeatureType TYPE = builder.buildFeatureType();
		for (int i=0; i<geometries.size(); i++) {
			int id = ids.get(i);
			Geometry geom = geometries.get(i);
			SimpleFeature feature = SimpleFeatureBuilder.build( TYPE, new Object[]{geom, id}, null );
			//				collection.add( feature );
			this.features.put(id, feature);
			this.geometries.put(id, geom);
		}
		return true;
	}

	/**
	 * Write out all the features stored in this Shapefile to a shapefile.
	 * @param outputFileName The name of the shapefile to write to.
	 * @return True if the operation was a success, false otherwise.
	 */
	public boolean writeShapefile(String outputFileName) {
		if (!this.checkFeatures()) {
			System.err.println("Shapefile.writeShapefile() error: no features to be written");
			return false;
		}
		File outFile = new File(outputFileName);
		// Check the path to the new file exists
		File parent = new File(outFile.getParent()); 
		if ( !parent.isDirectory() ) {
			parent.mkdirs();
		}
		
		// Check the output file can be written to
		//		if (!outFile.canWrite()) {
		//			System.err.println("Shapefile.writeShapefile() error: cannot write to the shapefile: "+outputFileName);
		//			return false;
		//		}		
		// Create a feature collection to write the features out to
		FeatureCollection<SimpleFeatureType, SimpleFeature> outFeatures = FeatureCollections.newCollection();
		for (SimpleFeature f:this.features.values()) {
			outFeatures.add(f);
		}
		try {
			// Don't really get the rest, copied from the GeoTools tutorial.
			DataStoreFactorySpi factory = new ShapefileDataStoreFactory();
			Map<String, Serializable> create = new HashMap<String, Serializable>();
			create.put("url", outFile.toURI().toURL());
			create.put("create spatial index", Boolean.TRUE);
			ShapefileDataStore newDataStore = (ShapefileDataStore) factory.createNewDataStore(create);
			newDataStore.createSchema(outFeatures.getSchema());
			if (this.crs != null) {
				newDataStore.forceSchemaCRS(this.crs);
			}
			Transaction transaction = new DefaultTransaction("create");
			String typeName = newDataStore.getTypeNames()[0];
			FeatureStore<SimpleFeatureType, SimpleFeature> featureStore;
			featureStore = (FeatureStore<SimpleFeatureType, SimpleFeature>) newDataStore.getFeatureSource(typeName);
			featureStore.setTransaction(transaction);
			try {
				featureStore.addFeatures(outFeatures);
				transaction.commit();
			} catch (Exception problem) {
				System.out.println("Shapefile.writeShapefile() caught a problem trying to write: "+problem.toString());
				problem.printStackTrace();
				transaction.rollback();
				return false;
			} finally {
				transaction.close();
			}

		} catch (IOException e) {
			System.err.println("Shapefile.writeShapefile() caught an IOException trying to write shapefile.");
			e.printStackTrace();
			return false;
		}
		return true;
	}

	/**
	 * Clears all objects from this Shapefile, useful for when the shape file is going to be re-read.
	 */
	private void clear() {
		this.geometries.clear();
		this.features.clear();
	}

	/**
	 * Resets this <code>Shapefile</code> by removing all existing features and re-reading the original
	 * shapefile.
	 */
	public void reset() {
		this.clear();
		this.readShapefile();
	}

	/**
	 * Return the ID's of all the features in this shapfile. This can be used to iterate over all
	 * <code>SimpleFeature</code>s or all <code>Geometry</code>s.
	 * @return The ID of every feature currently in this <code>Shapefile</code>.
	 */
	public Set<Integer> getFeatureIDs() {
		return this.features.keySet();
	}

	/**
	 * Get the feature with the associated ID.
	 * @param id
	 * @return The Feature or null if no feature is found.
	 */
	public SimpleFeature getFeature(int id) {
		if (!this.features.containsKey(id)) {
			System.out.println("Shapefile.getFeature() error, no feature found with ID: "+id);
		}
		return this.features.get(id);
	}

	public List<SimpleFeature> getFeatures() {
		List<SimpleFeature> out = new ArrayList<SimpleFeature>(this.features.values().size());
		for (SimpleFeature f:this.features.values()) {
			out.add(f);
		}
		return out;
	}

	/**
	 * Get the geometry of the object with the associated ID.
	 * @param id
	 * @return The Geometry or null if no feature is found.
	 */
	public Geometry getGeometry(int id) {
		if (!this.geometries.containsKey(id)) {
			System.out.println("Shapefile.getFeature() error, no feature found with ID: "+id);
		}
		return this.geometries.get(id);
	}

	/** 
	 * Get all the geometries of the objects.
	 * @return A list of geometries or null if none are found.
	 */
	public List<Geometry> getGeometries() {
		List<Geometry> l = new ArrayList<Geometry>(this.geometries.values().size());
		for (Geometry g:this.geometries.values()) {
			l.add(g);
		}
		return l;
	}

	/**
	 * Get the name of this <code>Shapefile</code>.
	 */
	public String getName() {
		return this.file.getName();
	}

	/**
	 * Create a clone of this Shapefile.
	 */
	@SuppressWarnings("unchecked")
	@Override
	protected Object clone() throws CloneNotSupportedException {
		Shapefile s = new Shapefile();
		s.geometries = (Hashtable<Integer, Geometry>) this.geometries.clone();
		s.features = (Hashtable<Integer, SimpleFeature>) this.features.clone();
		s.idColumn = this.idColumn;
		s.file = new File(this.file.toURI());
		return s;
	}


}
