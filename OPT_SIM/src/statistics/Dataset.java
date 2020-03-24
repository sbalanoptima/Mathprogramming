/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package statistics;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.Scanner;

/**
 *
 * @author MBoon
 */
public class Dataset {

  private String[] varNames; // variable names
  private double[][] values; // values


  
  /**
   * Imports a dataset from a tab delimited text file. The first line should contain the variable name(s)
   * @param inputFile the input file
   */
  
  public Dataset(File inputFile) throws IOException {
      Scanner fileScanner = new Scanner(inputFile);
      // read the first line
      try {
        String line = fileScanner.nextLine();
        Scanner stringScanner = new Scanner(line);
        ArrayList<String> colNames = new ArrayList<String>();
        while (stringScanner.hasNext()) {
            colNames.add(stringScanner.next());
        }
        if (colNames.size() == 0) {
            throw new IOException("Error: the first line does not have any entries.");
        }
        this.varNames = new String[colNames.size()];
        for (int i = 0; i < varNames.length; i++) {
            varNames[i] = colNames.get(i);
        }
      } catch (NoSuchElementException e) {
          throw new IOException("Error: the selected file appears to be empty.");
      }
      ArrayList<double[]> v = new ArrayList<double[]>();
      while (fileScanner.hasNext()) {
          String line = fileScanner.nextLine();
          Scanner stringScanner = new Scanner(line);
          double[] lineValues = new double[varNames.length];
          int i = 0;
          while (stringScanner.hasNext() && i < varNames.length) {
              lineValues[i] = Double.parseDouble(stringScanner.next());
              i++;
          }
          v.add(lineValues);
      }
      values = new double[v.size()][varNames.length];
      for (int i = 0; i < values.length; i++) {
          values[i] = (double[])v.get(i);
      }
      fileScanner.close();
  }
  

  public Dataset(String[] variableNames, double[][] values) {
      this.varNames = new String[variableNames.length];
      this.values = new double[values.length][values[0].length];
      for (int i = 0; i < variableNames.length; i++) {
          this.varNames[i] = variableNames[i];
          for (int j = 0; j < values.length; j++) {
              this.values[j][i] = values[j][i];
          }
      }
  }

  public double[][] getValues() {
      return values;
  }
  
  public String[] getVariableNames() {
      return varNames;
  }
  
  public void export(File fileName) throws IOException {
          PrintWriter pw = new PrintWriter(new FileWriter(fileName));
          for (int i = 0; i < varNames.length; i++) {
              pw.print(varNames[i] + "\t");
          }
          pw.println();
          for (int j = 0; j < values.length; j++) {
              for (int i = 0; i < values[j].length; i++) {
                  pw.print(values[j][i]);
                  if (i < values[j].length - 1) { // not yet end of line -> next column
                      pw.print("\t"); // next column
                  } else {
                      pw.println();  // next line
                  }
              }
          }
          pw.flush();
          pw.close();
  }
  
  public String toString() {
      String s = "";
      for (int i = 0; i < varNames.length; i++) {
          s += varNames[i];
          if (i < varNames.length - 1) {
              s += "\t";
          } else {
              s += "\n";
          }
      }
      for (int j = 0; j < values.length; j++) {
          for (int i = 0; i < varNames.length; i++) {
              s += values[j][i];
              if (i < varNames.length - 1) {
                  s += "\t";
              } else {
                  s += "\n";
              }
          }
      }
      return s;
  }
  
  public static void main(String[] arg) {
      String[] vars = {"Stud.Id.", "Grade"};
      double[][] values = {
          {345346, 8.5},
          {399596, 9.1},
          {435231, 6.5},
          {234565, 8.2}
      };
      Dataset data = new Dataset(vars, values);
      try {
        File testFile = new File("D:/temp/grades.txt");
        data.export(testFile);
        Dataset data2 = new Dataset(testFile);
        System.out.println(data2);
      } catch (IOException e) {
          e.printStackTrace();
      }
  }
}
