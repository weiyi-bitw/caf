package graphicTest;

import javax.swing.*;

import util.WindowUtilities;

import java.awt.*;
import java.awt.geom.*;

public class ShapeExample extends JPanel {
	private Ellipse2D.Double circle =
		    new Ellipse2D.Double(10, 10, 350, 350);
	private Rectangle2D.Double square =
		    new Rectangle2D.Double(10, 10, 350, 350);

	public void paintComponent(Graphics g) {
		    clear(g);
		    Graphics2D g2d = (Graphics2D)g;
		    g2d.fill(circle);
		    g2d.draw(square);
	}

	// super.paintComponent clears offscreen pixmap,
	// since we're using double buffering by default.

	protected void clear(Graphics g) {
		super.paintComponent(g);
	}

	protected Ellipse2D.Double getCircle() {
	    return(circle);
	}

	public static void main(String[] args) {
	    WindowUtilities.openInJFrame(new ShapeExample(), 380, 400);
	}	
	
}
