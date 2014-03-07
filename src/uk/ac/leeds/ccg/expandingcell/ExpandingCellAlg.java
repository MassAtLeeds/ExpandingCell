package uk.ac.leeds.ccg.expandingcell;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

import javax.swing.JTextArea;

import org.opengis.feature.simple.SimpleFeature;

import uk.ac.leeds.ccg.expandingcell.Shapefile;
import uk.ac.leeds.mass.fmf.fit_statistics.GOF_TEST;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.Polygon;

/*
 * Implements Andy's expanding cell error measuring algorithm.
 * 
 */
public class ExpandingCellAlg implements Runnable {

    private List<SimpleFeature> expectedFeatures; // Feature lists
    private List<SimpleFeature> predictedFeatures;
    private List<SimpleFeature> allFeatures;		// Combination of both lists
    private Cell box; // The bounding box
    // Gathers results from each iteration. 2D list, stores number of results at each resolution (grid shifted)
    private List<List<Results>> results;
    private int increment;
    private int upperLimit; // Max number of cells per row
    private int outputTime; // The point at which to output a grid of results (num cells per row)
    //	private GeometryFactory geomFac;
    private JTextArea textArea;	// For sending output to the gui. Will use System.out if this is null
    private File expShapeF, preShapeF, outShapeF, outCSVF; // In/Out files

    /**
     * Create a new algorithm object, ready to run.
     *
     * @param expectedShapefile The shape file with the original data (expected point locations)
     * @param predictedShapefile The shape file with the model generated data (predicted point locations)
     * @param outputShapeFile Where to write out a shapefile of the results to.
     * @param outputCSVFile Where to write out all the results of each grid size as a csv file.
     * @param upperLimit Maximum number of cells per row
     * @param increment How many new cells per row to add after each increment (can be used to speed up algorithm).
     * @param outputTime The point at which to output the results to a shapefile (iteration number). I.e. run
     * the algorithm, graph the errors, find an interesting point on the graph and then run the algorith again,
     * this time outputting a shapefile of the results from that point.
     * @param textArea An optional TextArea to send results to. This can be used in conjunction with a GUI, if
     * it is null then output is sent to console.
     */
    public ExpandingCellAlg(File expectedShapefile, File predictedShapefile,
            File outputShapeFile,
            File outputCSVFile, int upperLimit, int increment, int outputTime,
            JTextArea textArea) {
        this.expShapeF = expectedShapefile;
        this.preShapeF = predictedShapefile;
        this.outShapeF = outputShapeFile;
        this.outCSVF = outputCSVFile;
        this.upperLimit = upperLimit;
        this.increment = increment;
        this.outputTime = outputTime;
        this.textArea = textArea;
    }

    /**
     * Run the algorithm.
     */
    public void run() {
        this.output("Running model with upper limit: " + this.upperLimit + " and increment: " + increment + "\n");
        expectedFeatures = new ArrayList<SimpleFeature>();
        predictedFeatures = new ArrayList<SimpleFeature>();
        allFeatures = new ArrayList<SimpleFeature>();

        // Create Shapefile objects and read in the features
        Shapefile expShape = new Shapefile(expShapeF);
        expShape.readShapefile();
        Shapefile preShape = new Shapefile(preShapeF);
        preShape.readShapefile();
        expectedFeatures = expShape.getFeatures();
        predictedFeatures = preShape.getFeatures();

        System.out.println("Found " + expectedFeatures.size() + " expected and "
                + predictedFeatures.size() + " predicted features.");

        // Combine them and calculate the bounding box
        allFeatures.addAll(expectedFeatures);
        allFeatures.addAll(predictedFeatures);
        this.box = calcBoundingBoxCoords(allFeatures);

        results = new ArrayList<List<Results>>(); // Store result of each cell size here.
        calculateErrors(results); // *MAIN ALGORITHM* populate results array

        // Print results:
        String res = "*Results (in csv format)* (also written to " + this.outCSVF.getName() + ")\n";
        res += Results.getCSVHeader() + "\n";
        for (int i = 0; i < results.size(); i++) {
            for (int j = 0; j < results.get(i).size(); j++) {
                res += results.get(i).get(j).toString() + "\n";
            } // for j
        } // for i

        this.output(res);
        // Write to text file:
        this.output("Writing results to csv file...");
        writeToCSVFile();
        this.output("...done\n");

        if (outShapeF == null) {
            System.out.println("Output shapfile is null so not writing any shapefiles");
        }
        else {
            // If outputTime is -1 then print out the last set of results (smallest grid size)
            // TODO Error here, if increment is > 1 then results.get(outputTime) will throw nullpointer. Need to
            // go through results array to work out which index holds results for outputTime..
            writeFeaturesToShapefile(outShapeF, results.get(outputTime == -1 ? results.size() - 1 : outputTime));
        }

    }

    private void calculateErrors(List<List<Results>> theResults) {
        // Loop through each cell size and calculate the difference in the number of expected/predicted points.
        Cell currCell = null;			// The cell we're looking for points within
        Cell cellTemplate = null;		// The template for all cells of this size
        int currNumCellsRow = 1;		// THe number of cells in a row which will fit into the bounding box
        int totCells = 0;				// Total number of cells, only used for outputting inf
        int counter = 0;				// Counter for while loop (not important, for output more than anything else)

        while (currNumCellsRow < upperLimit) { 			// Iterate over each resolution
            //			// Create the new cell template which we'll use to search for points within
            //			cellTemplate = new Cell(box.minX, box.minY,
            //					(box.width()/currNumCellsRow) + box.minX,
            //					(box.height()/currNumCellsRow) + box.minY );
            // Cell template slightly outside sim area to allow for grid shifting
            cellTemplate = new Cell(box.minX, box.minY,
                    (box.width() / currNumCellsRow) + box.minX,
                    (box.height() / currNumCellsRow) + box.minY);
            totCells = (int) Math.pow(currNumCellsRow, 2);
            this.output("Iteration " + counter + ": cells are " + (1.0 / currNumCellsRow) * 100 + " percent of total size, " + totCells + " cells will be used.\n");
            theResults.add(new ArrayList<Results>()); // Create new list to hold results for this resolution
			/* Shift the entire grid four times (N, E, S, W). This is done by creating a grid slightly larger than
             * required (the cell template is shifted slightly SW and extra cells are created at end of column and row).*/
            for (int shift = 0; shift < 5; shift++) { 	//
                // Create the new Results object for this shift and add it to the list of all results
                Results currentResults = new Results(counter, totCells);
                theResults.get(theResults.size() - 1).add(currentResults);
                // Create a cell template which is shifted in the correct direction
                Cell currentCellTemplate = calcShiftedCellTempate(shift, cellTemplate);
                // Loop over every cell in the grid for this cell size

                // Have to create extra cells if shifting the grid, otherwise some points will be outsidegrid
                int numCellsToCreate = (shift == 0 ? currNumCellsRow : currNumCellsRow + 1);
                double minX, minY, maxX, maxY; // min/max x,y for currCell
                for (int i = 0; i < numCellsToCreate; i++) {
                    for (int j = 0; j < numCellsToCreate; j++) {
                        minX = currentCellTemplate.minX + (i * (currentCellTemplate.width())); 		// minX
                        minY = currentCellTemplate.minY + (j * (currentCellTemplate.height()));		// minY
                        maxX = minX + cellTemplate.width();				// maxX
                        maxY = minY + cellTemplate.height();				// maxY
                        currCell = new Cell(minX, minY, maxX, maxY);
                        //					System.out.println("\tcell coords:"+currCell.toString());
                        // Look for expected and predicted features, saving the results
                        //  (very inefficient, loops through every point for each cell!)
                        double exp = 0;
                        double pre = 0;
                        Coordinate coord;
                        for (SimpleFeature f : expectedFeatures) {
                            coord = ((Geometry) f.getDefaultGeometry()).getCoordinate();
                            if (within(coord, currCell)) {
                                currCell.expected += 1;
                                exp++;
                            }
                        }
                        for (SimpleFeature f : predictedFeatures) {
                            coord = ((Geometry) f.getDefaultGeometry()).getCoordinate();
                            if (within(coord, currCell)) {
                                currCell.predicted += 1;
                                pre++;
                            }
                        }
                        // Add this cell to the current results object
                        currentResults.addCell(currCell);
                    } // for j
                } // for i
                // Now have total expected and predicted numbers can calculate relative cell errors
                for (Cell c : currentResults.getCells()) {
                    double prePct = 100 * ((double) c.predicted / (double) currentResults.getTotalPredicted());
                    double expPct = 100 * ((double) c.expected / (double) currentResults.getTotalExpected());
                    c.PctExpected = expPct;
                    c.PctPredicted = prePct;
                    c.relPctError = prePct - expPct;
                    c.pctError = Math.abs(prePct - expPct);
                }
                //			System.out.println("\t Total points found: predicted: "+results.get(counter).getTotalPredicted()+
                //			", expected: "+results.get(counter).getTotalExpected()+". Total error: "+results.get(counter).calcError());
                //				this.output(". "+currentResults.toString()+"\n");
            } // for shift
            currNumCellsRow += this.increment;
            counter++;
        } // while
    }

    /**
     * Takes a cellTemplate (the south-western cell in the corner of the bounding box) and returns a "shifted"
     * Cell.
     *
     * @param shift An integer: 0 - no shift, 1 - north, 2 - east, 3 - south, 4 - west
     * @param cellTemplate The original cell template (lower-left corner of bounding box)
     * @return Shifted cell
     */
    private Cell calcShiftedCellTempate(int shift, Cell cellTemplate) {
        // First need to move template cell south-west a bit (half cell size height)
        Cell theTemplate = new Cell(
                cellTemplate.minX - ((double) cellTemplate.width() / 2.0),
                cellTemplate.minY - ((double) cellTemplate.height() / 2.0),
                cellTemplate.maxX - ((double) cellTemplate.width() / 2.0),
                cellTemplate.maxY - ((double) cellTemplate.height() / 2.0));

        // Amount to shift cell grid each time (quarter or the cell size)
        double shiftXAmount = ((double) cellTemplate.width()) / 4.0;
        double shiftYAmount = ((double) cellTemplate.height()) / 4.0;

        // Now shift cell depending on the direction
        if (shift == 0) { 			// No shift
            return cellTemplate;
        }
        else if (shift == 1) { // Shift quarter of a cell north
            return new Cell(
                    theTemplate.minX,
                    theTemplate.minY + shiftYAmount,
                    theTemplate.maxX,
                    theTemplate.maxY + shiftYAmount);
        }
        else if (shift == 2) { // Shift quarter of a cell east
            return new Cell(
                    theTemplate.minX + shiftXAmount,
                    theTemplate.minY,
                    theTemplate.maxX + shiftXAmount,
                    theTemplate.maxY);
        }
        else if (shift == 3) { // Shift quarter of a cell south
            return new Cell(
                    theTemplate.minX,
                    theTemplate.minY - shiftYAmount,
                    theTemplate.maxX,
                    theTemplate.maxY - shiftYAmount);
        }
        else if (shift == 4) { // Shift quarter of a cell west
            return new Cell(
                    theTemplate.minX - shiftXAmount,
                    theTemplate.minY,
                    theTemplate.maxX - shiftXAmount,
                    theTemplate.maxY);
        }
        else {
            this.output("ExpandingCell.calcShiftedCellTemplate error, unrecognised shift value: " + shift);
            return null;
        }
    }

    private void output(String in) {
        if (this.textArea == null) {
            System.out.print(in);
        }
        else {
            textArea.append(in);
        }
    }

    private Cell calcBoundingBoxCoords(List<SimpleFeature> features) {
        double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE, maxX = Double.NEGATIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY;
        double x, y;
        int i = 0;
        for (SimpleFeature f : features) {
            x = ((Geometry) f.getDefaultGeometry()).getCoordinate().x;
            y = ((Geometry) f.getDefaultGeometry()).getCoordinate().y;
            if (x < minX) {
                minX = x;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
            i++;
        }
        //		// Now make slightly larger (0.11%) to make sure we get every point
        //		Cell box = new Cell(minX-(minX*0.001), minY-(minY*0.001), maxX+(maxX*0.001), maxY+maxY*0.001);
        Cell box = new Cell(minX, minY, maxX, maxY);
        this.output("Found LL/UR bounding box coords: (" + box.minX + "," + box.minY + ") / (" + box.maxX + "," + box.maxY + ")\n");
        return box;

    }

    // Is Coordinate c within Cell cell
    private boolean within(Coordinate c, Cell cell) {
        // XXXX if point lies on boundary it will be counted twice
        return (c.x <= cell.maxX) && (c.x >= cell.minX) && (c.y <= cell.maxY) && (c.y >= cell.minY);
    }

    private void writeFeaturesToShapefile(File outFile, List<Results> theResults) {
        this.output("Writing features to shapefile...");
        // Run through each shifted grid
        for (int shift = 0; shift < theResults.size(); shift++) {
            // Add the shift number to the shapefile (input shapefile doesn't include it)
            File theOutFile = new File(outFile.getAbsolutePath().substring(0, outFile.getAbsolutePath().length() - 3) + "-" + shift + ".shp");
            Shapefile outCell = new Shapefile(theOutFile);
            List<Integer> ids = new ArrayList<Integer>();	// IDs of each cell
            List<Geometry> geometries = new ArrayList<Geometry>();	// Geometries of each cell
            Map<Integer, Double> cellSizes = new Hashtable<Integer, Double>();	// The size of each cell
            Map<Integer, Double> absErrors = new Hashtable<Integer, Double>();	// The absolute error of each cell
            Map<Integer, Double> relPctErrors = new Hashtable<Integer, Double>();	// The relative % of each cell
            Map<Integer, Double> absPctErrors = new Hashtable<Integer, Double>();	// The absoluve % error of each cell
            Map<Integer, Integer> numExpected = new Hashtable<Integer, Integer>();	// The number of predicted points
            Map<Integer, Integer> numPredicted = new Hashtable<Integer, Integer>();	// The number of expected points
            Map<Integer, Double> pctExpected = new Hashtable<Integer, Double>();	// The % predicted points in cell
            Map<Integer, Double> pctPredicted = new Hashtable<Integer, Double>();	// The % expected points in cell
            GeometryFactory cellGeomFac = new GeometryFactory();

            // Run through each cell, writing it's geometry and error to the datastore
            Results currentResults = theResults.get(shift);
            for (int i = 0; i < currentResults.getCells().size(); i++) {
                Cell cell = currentResults.getCells().get(i);
                Coordinate[] ringCoords = new Coordinate[]{
                    new Coordinate(cell.minX, cell.minY), // LL
                    new Coordinate(cell.minX, cell.maxY), // UL
                    new Coordinate(cell.maxX, cell.maxY), // UR
                    new Coordinate(cell.maxX, cell.minY), // LR
                    new Coordinate(cell.minX, cell.minY), // LL
                };
                LinearRing shape = cellGeomFac.createLinearRing(ringCoords);
                Polygon geometry = cellGeomFac.createPolygon(shape, new LinearRing[0]);
                geometries.add(geometry);
                ids.add(i);
                cellSizes.put(i, cell.size);
                numExpected.put(i, cell.expected);
                numPredicted.put(i, cell.predicted);
                pctExpected.put(i, cell.PctExpected);
                pctPredicted.put(i, cell.PctPredicted);
                absErrors.put(i, cell.absError);
                relPctErrors.put(i, cell.relPctError);
                absPctErrors.put(i, cell.pctError);
            } // for cells

            outCell.createFeatures(geometries, ids, "ID", null);
            outCell.addAttribute("CellSize", Double.class, cellSizes);
            outCell.addAttribute("Predicted", Integer.class, numPredicted);
            outCell.addAttribute("Expected", Integer.class, numExpected);
            outCell.addAttribute("AbsError", Double.class, absErrors);
            outCell.addAttribute("PctPred", Double.class, pctPredicted);
            outCell.addAttribute("PctExpec", Double.class, pctExpected);
            outCell.addAttribute("AbsPctErr", Double.class, absPctErrors);
            outCell.addAttribute("RelPctErr", Double.class, relPctErrors);

            outCell.writeShapefile(theOutFile.getAbsolutePath());
        } // for shift
        this.output("done.\n");
    }

    private void writeToCSVFile() {
        try {
            // Check directorys don't need to be created
            File parent = new File(this.outCSVF.getParent());
            if (!parent.isDirectory()) {
                parent.mkdirs();
            }
            BufferedWriter w = new BufferedWriter(new FileWriter(this.outCSVF));
            w.write(Results.getCSVHeader());

            for (int i = 0; i < results.size(); i++) {
                String resultsStr = new String();
                for (Results r : results.get(i)) {
                    resultsStr += r.toString() + "\n";
                }
                w.write(resultsStr);
            }
            w.close();

        }
        catch (IOException e) {
            output("ExpandingCellAlg.writeToCSVFile IOException:");
            e.printStackTrace();
        }
    }

    /* Provided for convenience, should use GUI to start algorithm */
    public static void main(String[] args) {

        int outputTime = 100;

        if (args.length == 4) {
            System.out.println("Using command line arguments to start algorithm:"
                    + "\n\tInput 1: " + args[0]
                    + "\n\tInput 2: " + args[1]
                    + "\n\tShapefile out: " + args[2]
                    + "\n\tCSV out: " + args[3]);

            // Run on twitter data
            ExpandingCellAlg a = new ExpandingCellAlg(
                    new File(args[0]),
                    new File(args[1]),
                    new File(args[2]), // out shapefile
                    new File(args[3]), // out csv
                    outputTime + 2, // upper limit
                    10, // increment,
                    -1, // output time (iteration number)
                    null // text area
                    );
            a.run();

        }
        else {
            System.out.println("There haven't been four command-line arguments specified "
                    + "so I'm running the algorithm on hard-coded data");

            // Run on twitter data
            ExpandingCellAlg a = new ExpandingCellAlg(
                    new File("/Volumes/RESDATA/research-data/projects/twitter/1/temporal_hotspots/gisruk/afternoons/afternoons.shp"),
                    new File("/Volumes/RESDATA/research-data/projects/twitter/1/temporal_hotspots/gisruk/evenings/evenings.shp"),
                    new File("/Volumes/RESDATA/research-data/projects/twitter/1/temporal_hotspots/gisruk/expanding_cell/ec-aft-eve-" + outputTime + ".shp"), // out shapefile
                    new File("/Volumes/RESDATA/research-data/projects/twitter/1/temporal_hotspots/gisruk/expanding_cell/ec-aft-eve-" + outputTime + ".csv"), // out csv
                    outputTime + 2, // upper limit
                    10, // increment,
                    -1, // output time (iteration number)
                    null // text area
                    );
            a.run();

        }

    }

}

class Cell {

    public double minX;
    public double minY;
    public double maxX;
    public double maxY;
    public double size; 		// The size (in square units) of the cell
    public int expected;	// The number of predicted and expected points within this cell
    public int predicted;
    public double PctExpected; 	// Percentage of expected points in this cell
    public double PctPredicted;	// Percentage of predicted points in this cell
    public double absError;	// Absolute difference in number of predicted - expected points
    public double pctError;	// Absolute percentage error
    public double relPctError;	// Relative percentage difference (PctExpected - PctExpected)

    public Cell(double minX, double minY, double maxX, double maxY) {
        this.minX = minX;
        this.minY = minY;
        this.maxX = maxX;
        this.maxY = maxY;
        this.size = (this.maxX - this.minX) * (this.maxY - this.minY);
    }

    public double width() {
        return this.maxX - this.minX;
    }

    public double height() {
        return this.maxY - this.minY;
    }

    public String toString() {
        return (this.minX + "," + this.minY + "," + this.maxX + "," + this.maxY);
    }

}

/**
 * Used to store results from an iteration and calculate the total error.
 * @author Nick Malleson
 */
class Results {

    public int iterationNum; 		// The iteration number;
    public int cellSize;			// Size of cells (as percent of total bounding box size)
    public double cellSizeAbs = -1;	// Size of cells (as squre units)
    private ArrayList<Cell> cells;	// All the cells
    public double error;			// The total error of this iteration
    private int totalPredicted = 0;
    private int totalExpected = 0;

    public Results(int iterationNum, int cellSize) {
        this.cells = new ArrayList<Cell>();
        this.iterationNum = iterationNum;
        this.cellSize = cellSize;
    }

    /**
     * Calculate SRMSE error.
     */
    public double calcSRMSE() {
        // Need to produce a double[][] matrix of expected and predicted cell values to pass to GOF test
        double[][] calib = new double[1][cells.size()];
        double[][] test = new double[1][cells.size()];
        for (int i = 0; i < cells.size(); i++) {
            calib[0][i] = (double) cells.get(i).expected;
            test[0][i] = (double) cells.get(i).predicted;
        }
        // Normalise the matricies so they both sum to 1
        GOF_TEST.normalise(calib);
        GOF_TEST.normalise(test);

        return GOF_TEST.SRMSE.getTest().test(calib, test);
    }

    /**
     * Calculate R2 (r-squared) error
     * @return
     */
    public double calcR2Error() {

        // Need to produce a double[][] matrix of expected and predicted cell values to pass to GOF test
        double[][] calib = new double[1][cells.size()];
        double[][] test = new double[1][cells.size()];
        for (int i = 0; i < cells.size(); i++) {
            calib[0][i] = (double) cells.get(i).expected;
            test[0][i] = (double) cells.get(i).predicted;
        }
        // Normalise the matricies so they both sum to 1
        GOF_TEST.normalise(calib);
        GOF_TEST.normalise(test);
        return GOF_TEST.R2.getTest().test(calib, test);
    }

    /**
     * Calculate the absolute error ( |predicted - expected| / predicted )
     * @return
     */
    public double calcTotalAbsError() {
        double totalError = 0;
        for (Cell c : cells) {
            if (!(c.predicted == 0 && c.expected == 0)) { // Cell must have a data point in it
                if (c.predicted == 0) {
                    totalError += Math.abs(c.expected - c.predicted); // Don't divid by zero
                }
                else {
                    totalError += Math.abs(c.expected - c.predicted) / c.predicted;
                }
                //				totalError += Math.sqrt(Math.pow((c.expected-c.predicted),2));

            }

        }
        return (totalError);
    }

    public void addCell(Cell cell) {
        cells.add(cell);
        if (this.cellSizeAbs == -1) { // If the cell size hasn't been set, set it
            this.cellSizeAbs = cell.size;
        }
        this.totalPredicted += cell.predicted;
        this.totalExpected += cell.expected;
        // Tell this cell what it's absolute error is (can be included in a shapefile)
        cell.absError = Math.abs(cell.predicted - cell.expected);
        //		cell.error = ( cell.predicted==0 ? // Allow for divide by zero
        //				(Math.abs(cell.expected-cell.predicted)) :
        //				(Math.abs(cell.expected-cell.predicted)/cell.predicted)
        //				);
    }

    public ArrayList<Cell> getCells() {
        return cells;
    }

    public int getTotalPredicted() {
        return this.totalPredicted;
    }

    public int getTotalExpected() {
        return this.totalExpected;
    }

    /**
     * Get the iteration number, cell size and error (in csv format).
     */
    public String toString() {
        DecimalFormat format = (DecimalFormat) NumberFormat.getInstance();
        format.applyPattern("#######.#################");
        // NOTE: convert the cellSize (unknown units) to square km by multiplying by 20970.55793
        // (found this number by measuring resulting area cell in a GIS)
        return iterationNum + ", " + cellSize + ", " + this.cellSizeAbs + ", " + (this.cellSizeAbs * 20970.55793) + ","
                + format.format(calcTotalAbsError()) + "," + format.format(calcR2Error()) + "," + format.format(calcSRMSE());
    }

    public static String getCSVHeader() {
        return "Iteration, TotCells, CellSize, CellSize_km2, TotalAbsoluteError, R2Error, SRMSE\n";
    }

}
