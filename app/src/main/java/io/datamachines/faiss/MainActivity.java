package io.datamachines.faiss;

import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.SeekBar;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity {

    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }

    EditText editor;
    SeekBar bar;
    Button StartTest;
    TextView tv;
    TextView tv2;
    TextView tv3;
    TextView tv4;
    TextView tv5;
    TextView tv6;
    String data;
    int count = 1;
    int currentValue;
    void StartTesting(){
        StartTest.setEnabled(false);
        tv.setText("Testing for now - Build time");

        tv2.setText("Testing for now - Train time");

        tv3.setText("Testing for now - dB Time");

        tv4.setText("Testing for now - Search Time");

        tv5.setText("Testing for now - Total Size");

        tv6.setText("Testing for now - Counter");
        currentValue = bar.getProgress()*10000;

        editor.setText("Testing with size ="+currentValue);
       new Thread(new Runnable()
       {
           @Override
           public void run()
           {
               try
               {


                   data = stringFromJNI(currentValue);
                   runOnUiThread(new Runnable()
                       {
                           @Override public void run()
                           {
                               UpdateGUI();
                           }
                       });

               }
               catch (Exception e)
               {

               }
           }
       }).start();

    }
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        // Example of a call to a native method
        editor = findViewById(R.id.editText);
        tv = findViewById(R.id.sample_text);
        tv.setText("Testing for now - Build time");

        tv2 = findViewById(R.id.sample_text2);
        tv2.setText("Testing for now - Train time");

        tv3 = findViewById(R.id.sample_text3);
        tv3.setText("Testing for now - dB Time");

        tv4 = findViewById(R.id.sample_text4);
        tv4.setText("Testing for now - Search Time");

        tv5 = findViewById(R.id.sample_text5);
        tv5.setText("Testing for now - Total Size");
        tv6 = findViewById(R.id.sample_text6);

        tv6.setText("Testing for now - Counter");
        TextView tv5;

        bar = findViewById(R.id.seekBar);


        StartTest = findViewById(R.id.button);
        StartTest.setOnClickListener(new View.OnClickListener() {
                              public void onClick(View v) {
                                  StartTesting();
                              }
                          });

        new Thread(new Runnable()
        {
            @Override
            public void run()
            {
                try
                {

                    while(true) {
                        Thread.sleep(2500);
                        runOnUiThread(new Runnable() {
                            @Override
                            public void run() {
                                UpdateCnt()     ;

                            }
                        });
                        count++;
                    }



                }
                catch (Exception e)
                {

                }
            }                                             
        }).start();

        bar.setOnSeekBarChangeListener(new SeekBar.OnSeekBarChangeListener() {
            @Override
            public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
                editor.setText("Testing size at "+progress*10000);
            }

            @Override
            public void onStartTrackingTouch(SeekBar seekBar) {

            }

            @Override
            public void onStopTrackingTouch(SeekBar seekBar) {

            }


        } );



























    }

    public void UpdateGUI()
    {
        String[] token = data.split(",");
        tv.setText(token[0]);
        tv2.setText(token[1]);
        tv3.setText(token[2]);
        tv4.setText(token[3]);
        tv5.setText(token[4]);

        StartTest.setEnabled(true);
    }

    public void UpdateCnt(){
         tv6.setText("Counter is at "+count);

    }

    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    public native String stringFromJNI(int a);
}
