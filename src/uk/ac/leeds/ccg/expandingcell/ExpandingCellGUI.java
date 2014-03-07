package uk.ac.leeds.ccg.expandingcell;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import java.io.File;


public class ExpandingCellGUI extends JPanel implements ActionListener {

	private static final long serialVersionUID = 1L;
	
	private JButton inExpB, inPreB, outShapeB, outCSVB, runB, defaultRunB;
	private File expShapeF, preShapeF, outShapeF, outCSVF;
	private JLabel expShapeLab, preShapeLab, outShapeLab, outCSVLab; // labels showing which files have been selected
	private JPanel buttonPanel, textPanel; // Separate panel for buttons 
    private JTextArea log;		// Reports from the console
    private JTextField uLimitT, incrementT, outputTimeT;	// Set the upper limit (cells per row) and how much to increment by
    private int uLimit = 100, increment = 1, outputTime = -1;
    private JFileChooser shapefc, savefc;	// For choosing shapefiles and directories (for saving) respectively
    private JScrollPane logScrollPane; 
 

    public ExpandingCellGUI() {
    	super (new BorderLayout());
    	this.setPreferredSize(new Dimension(600,800));
    	// Default values for filenames
//    	this.expShapeF = new File("data/burglaries-predicted_point.shp");
//    	this.preShapeF = new File("data/burglaries-expected_point.shp");
    	this.expShapeF = new File("data/observed.shp");
    	this.preShapeF = new File("data/simulated.shp");
    	this.outShapeF = new File("results/results.shp");
    	this.outCSVF = new File("results/results.csv");

        this.createTextFields();
        this.createFileChoosers();
        this.createButtons();

        //Add everything to the panel.
        add(buttonPanel, BorderLayout.PAGE_START);
        // (textPanel is actually added to the bottom of the buttonPanel - looks better this way)
        //add(textPanel, BorderLayout.LINE_END); 
        add(logScrollPane, BorderLayout.CENTER);
    }
    
    private void createTextFields() {
        //Create the log first, because the action listeners
        //need to refer to it.
        log = new JTextArea(5,20);
        log.setMargin(new Insets(5,5,5,5));
        log.setEditable(false);
        logScrollPane = new JScrollPane(log);
        
        // Create text fields for upperlimit and increment using default values
        uLimitT = new JTextField(String.valueOf(uLimit),10);
        uLimitT.addActionListener(this);
        uLimitT.setToolTipText("Terminate when this number of cells per row has been reached");
        incrementT = new JTextField(String.valueOf(increment),10);
        incrementT.addActionListener(this);
        incrementT.setToolTipText("Set the number of cells per row to add after each iteration");
        outputTimeT = new JTextField(String.valueOf(outputTime),10);
        outputTimeT.addActionListener(this);
        outputTimeT.setToolTipText("Output shapefile grid when this number of cells per row has been reached"); 
        textPanel = new JPanel(new GridLayout(3,2));
        // labels for the text fields
        JLabel uLimitLab = new JLabel("Upper Limit (cells per row)");
        uLimitLab.setToolTipText("Terminate when this number of cells per row has been reached");
        JLabel incrementLab = new JLabel("Increment");
        incrementLab.setToolTipText("Set the number of cells per row to add after each iteration");
        JLabel outputTimeLab = new JLabel("Output Grid Time");
        outputTimeLab.setToolTipText("Output shapefile grid when this number of cells per row has been reached");
        textPanel.add(uLimitLab);
        textPanel.add(uLimitT);
        textPanel.add(incrementLab);
        textPanel.add(incrementT);
        textPanel.add(outputTimeLab);
        textPanel.add(outputTimeT);
    }
    
    private void createFileChoosers() {
        //Create file choosers
        shapefc = new JFileChooser();
        shapefc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        shapefc.setCurrentDirectory(new File("."));
        shapefc.setFileFilter(new ShapeFilter()); // Custom implementation of FileFilter
        savefc = new JFileChooser();
        savefc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        savefc.setCurrentDirectory(new File("."));
    }
    
    // Create buttons and lables which show if a file has been selected.
    private void createButtons() {
        inExpB = new JButton("Open expected (real) shapefile");
        inExpB.addActionListener(this);
        inPreB = new JButton("Open predicted (model) shapefile");
        inPreB.addActionListener(this);
        outShapeB = new JButton("Choose location to save error results grid");
        outShapeB.addActionListener(this);
        outCSVB = new JButton("Choose location to save CSV results");
        outCSVB.addActionListener(this);
        runB = new JButton("Run!");
        runB .addActionListener(this);
        defaultRunB = new JButton("(Run using default file locations)");
        defaultRunB .addActionListener(this);
        
        // Create labels
        expShapeLab = new JLabel(expShapeF.getName());
        preShapeLab = new JLabel(preShapeF.getName());
        outShapeLab = new JLabel(outShapeF.getName());
        outCSVLab = new JLabel(outCSVF.getName());
        
        //For layout purposes, put the buttons in a separate panel
        buttonPanel = new JPanel(new GridLayout(6,2)); //use FlowLayout
        buttonPanel.add(inExpB);
        buttonPanel.add(expShapeLab);
        buttonPanel.add(inPreB);
        buttonPanel.add(preShapeLab);
        buttonPanel.add(outShapeB);
        buttonPanel.add(outShapeLab);
        buttonPanel.add(outCSVB);
        buttonPanel.add(outCSVLab);
        buttonPanel.add(runB);
        //buttonPanel.add(defaultRunB);
        buttonPanel.add(textPanel);
   }
    
    private void runAlg() {
    	// XXXX Do some checks
    	log.append("Running algorithm, this might take a while, you'll get no indication how long it will take!");
    	ExpandingCellAlg alg = new ExpandingCellAlg(this.expShapeF, this.preShapeF, this.outShapeF, this.outCSVF,
    			this.uLimit, this.increment, this.outputTime, this.log);
  
    		// XXXX progress bar doesn't work
//        	JProgressBar p = new JProgressBar();
//            p.setIndeterminate(true);
//            p.setVisible(true);
    	new Thread(alg).start(); // run the algorithm (in a new thread)
    	
    }
    private void runAlgDefault() {
    	this.expShapeF = new File("data/burglaries-predicted_point.shp");
    	this.preShapeF = new File("data/burglaries-expected_point.shp");
    	this.outShapeF = new File("results/results.shp");
    	this.outCSVF = new File("results/results.csv");
    	runAlg() ;    	
    }
 

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() {
    	//Create and set up the window.
    	JFrame frame = new JFrame("Expanding Cell Comparatorp");
    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    	//Add content to the window.
    	frame.add(new ExpandingCellGUI());

    	//Display the window.
    	frame.pack();
    	frame.setVisible(true);
    }
    
    

    // Uses buttons and the file dialogue to assign correct files to expShapeF, preShapeF, outShapeF, outCSVF;
    public void actionPerformed(ActionEvent e) {
    	try {
    		if (e.getSource() == inExpB) {
    			int returnVal =shapefc.showOpenDialog(ExpandingCellGUI.this);
    			if (returnVal == JFileChooser.APPROVE_OPTION) {
    				expShapeF = shapefc.getSelectedFile();
    				expShapeLab.setText(expShapeF.getName());
    				log.append("Selected expected (real data) results shapefile: " + expShapeF.getAbsolutePath() + ".\n");
    			} else {
    				log.append("Open command cancelled by user.\n");
    			}
    			log.setCaretPosition(log.getDocument().getLength());
    		} 
    		else if (e.getSource() == inPreB) {
    			int returnVal = shapefc.showOpenDialog(ExpandingCellGUI.this);
    			if (returnVal == JFileChooser.APPROVE_OPTION) {
    				preShapeF = shapefc.getSelectedFile();
    				preShapeLab.setText(preShapeF.getName());
    				log.append("Selected predicted (model) results shapefile: " + preShapeF.getAbsolutePath()+ ".\n");
    			} else {
    				log.append("Open command cancelled by user.\n");
    			}
    			log.setCaretPosition(log.getDocument().getLength());
    		}
    		else if (e.getSource() == outShapeB) {
    			int returnVal = savefc.showOpenDialog(ExpandingCellGUI.this);
    			if (returnVal == JFileChooser.APPROVE_OPTION) {
    				String dirName = savefc.getSelectedFile().getAbsolutePath();
    				outShapeF = new File(dirName+"/results.shp");
    				outShapeLab.setText(outShapeF.getPath());
    				log.append("Will save output shapefile (error grid) to: " + outShapeF.getAbsolutePath()+ ".\n");
    			} else {
    				log.append("Open command cancelled by user.\n");
    			}
    			log.setCaretPosition(log.getDocument().getLength());
    		} 
    		else if (e.getSource() == outCSVB) {
    			int returnVal = savefc.showOpenDialog(ExpandingCellGUI.this);
    			if (returnVal == JFileChooser.APPROVE_OPTION) {
    				String dirName = savefc.getSelectedFile().getAbsolutePath(); 
    				outCSVF = new File(dirName+"/results.csv");
    				outCSVLab.setText(outCSVF.getPath());
    				log.append("Will save csv results to " + outCSVF.getAbsolutePath()+ ".\n");
    			} else {
    				log.append("Open command cancelled by user.\n");
    			}
    			log.setCaretPosition(log.getDocument().getLength());
    		} 
    		else if (e.getSource() == runB) {
    			runAlg();
    		} 
    		else if (e.getSource() == defaultRunB) {
    			runAlgDefault();
    		}
    		else if (e.getSource() == uLimitT) {
    			uLimit = Integer.parseInt(uLimitT.getText());
    			log.append("Upper limit is "+uLimit+" cells per row.\n");
    		}
    		else if (e.getSource() == incrementT) {
    			increment = Integer.parseInt(incrementT.getText());
    			log.append("Increment is "+increment+" cell per row added each iteration.\n");
    		}
    		else if (e.getSource() == outputTimeT) {
    			outputTime = Integer.parseInt(outputTimeT.getText());
    			log.append("Output grid will be drawn when there are "+outputTime+" cells per row.\n");
    		}
    	}catch (NumberFormatException error) {
    		JOptionPane.showMessageDialog(this, "Please enter a positive number");
    	}
    }
    public static void main(String[] args) {
        //Schedule a job for the event dispatch thread:
        //creating and showing this application's GUI.
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                //Turn off metal's use of bold fonts
                UIManager.put("swing.boldMetal", Boolean.FALSE); 
                createAndShowGUI();
            }
        });
    } // main
} // class

// File filter to filter out shapefiles (*.shp)
class ShapeFilter extends javax.swing.filechooser.FileFilter {
    public boolean accept(File f) {
        return f.getName().toLowerCase().endsWith(".shp");
    }
    
    public String getDescription() {
        return "Shapefiles";
    }
}