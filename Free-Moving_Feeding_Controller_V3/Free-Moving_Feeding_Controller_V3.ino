#include <Wire.h>
#include <Servo.h>
#include <Adafruit_RGBLCDShield.h>
#include <utility/Adafruit_MCP23017.h>
#include <Adafruit_MotorShield.h>
#include "utility/Adafruit_MS_PWMServoDriver.h"

#define ServoPin 32
#define LimitSwitchPin_Z 34
#define LimitSwitchPin_X 33
#define LeverInsertPin 28
#define LeverPin 29
#define MouseDetectorPin1 30
#define MouseDetectorPin2 31
#define ACpowerTriggerPin 35

#define PiezoBuzzerPin 4

#define DIR_X 27
#define STEP_X 25
#define EN_X 26
#define CHOP_X 24
#define DIR_Z 37
#define STEP_Z 38
#define EN_Z 39
#define CHOP_Z 40

#define ServoPulse_Min 750
#define ServoPulse_Max 2250

int ServoAngle[2] = {1868, 1260};
float rpm_z = 120;
uint8_t microsteps_z = 16;
float rpm_x = 120;
uint8_t microsteps_x = 16;
float rpm_plate = 40;
int steps_z = 519;
int steps_pellet = 200;
int steps_x = 787;
int pellets = 24;
int n = 0;

unsigned long duration = 2400000;
unsigned long time_begin;

boolean UseLever = false;
boolean UseLowerPlate = true;
boolean TimeOut = false;

Adafruit_RGBLCDShield lcd = Adafruit_RGBLCDShield();
Servo Actuator;
Adafruit_MotorShield AFMS = Adafruit_MotorShield();
Adafruit_StepperMotor *StepperMotor_Plate_Upper = AFMS.getStepper(200, 2);
Adafruit_StepperMotor *StepperMotor_Plate_Lower = AFMS.getStepper(200, 1);

void setup() {
  pinMode(ACpowerTriggerPin, OUTPUT);
  AFMS.begin();
  StepperMotor_Plate_Lower->setSpeed(rpm_plate);
  StepperMotor_Plate_Lower->step(1, BACKWARD, MICROSTEP);
  StepperMotor_Plate_Lower->step(1, FORWARD, MICROSTEP);
  StepperMotor_Plate_Upper->setSpeed(rpm_plate);
  StepperMotor_Plate_Upper->step(1, BACKWARD, MICROSTEP);
  StepperMotor_Plate_Upper->step(1, FORWARD, MICROSTEP);
  pinMode(EN_X, OUTPUT);
  digitalWrite(EN_X, LOW);
  pinMode(CHOP_X, OUTPUT);
  digitalWrite(CHOP_X, HIGH);
  pinMode(DIR_X, OUTPUT);
  pinMode(STEP_X, OUTPUT);
  digitalWrite(STEP_X, HIGH);
  pinMode(EN_Z, OUTPUT);
  digitalWrite(EN_Z, LOW);
  pinMode(CHOP_Z, OUTPUT);
  digitalWrite(CHOP_Z, HIGH);
  pinMode(DIR_Z, OUTPUT);
  pinMode(STEP_Z, OUTPUT);
  digitalWrite(STEP_Z, HIGH);

  pinMode(LeverInsertPin, OUTPUT);
  digitalWrite(LeverInsertPin, HIGH);
  pinMode(LeverPin, INPUT);
  pinMode(MouseDetectorPin1, INPUT);
  digitalWrite(MouseDetectorPin1, HIGH);
  pinMode(MouseDetectorPin2, INPUT);
  digitalWrite(MouseDetectorPin2, HIGH);
  pinMode(PiezoBuzzerPin, OUTPUT);
  pinMode(LimitSwitchPin_Z, INPUT_PULLUP);
  pinMode(LimitSwitchPin_X, INPUT_PULLUP);

  Serial.begin(9600);

  // set up the LCD's number of columns and rows:
  lcd.begin(16, 2);
  uint8_t buttons;
  boolean Adjustment = false;
  Adjustment = user_input(String("Make some adjust"), String("ment? No"), String("ment? Yes"), Adjustment);
  digitalWrite(ACpowerTriggerPin, HIGH);

  Actuator.attach(ServoPin, ServoPulse_Min, ServoPulse_Max);
  if (Adjustment) {
    boolean ServoAdjust = false;
    ServoAdjust = user_input(String("Adjust Servo Ang"), String("le? No"), String("le? Yes"), ServoAdjust);
    if (ServoAdjust) {
      for (int i = 0; i < 2; i++) {
        Actuator.writeMicroseconds((ServoAngle[i]));
        delay(500);
        buttons = lcd.readButtons();
        while ((buttons & BUTTON_SELECT)) {
          buttons = lcd.readButtons();
        }
        lcd.clear();
        lcd.print("Adjust Angle ");
        lcd.print(i + 1);
        while (!(buttons & BUTTON_SELECT)) {
          buttons = lcd.readButtons();
          if ((buttons & BUTTON_RIGHT) || (buttons & BUTTON_LEFT)) {
            if (buttons & BUTTON_RIGHT) {
              ServoAngle[i]++;
              if (ServoAngle[i] > ServoPulse_Max) {
                ServoAngle[i] = ServoPulse_Max;
              }
            }
            if (buttons & BUTTON_LEFT) {
              ServoAngle[i]--;
              if (ServoAngle[i] < ServoPulse_Min) {
                ServoAngle[i] = ServoPulse_Min;
              }
            }
            Actuator.writeMicroseconds((ServoAngle[i]));
          }
        }
        lcd.clear();
        lcd.print("Angle ");
        lcd.print(i + 1);
        lcd.print(" = ");
        lcd.print(ServoAngle[i]);
        delay(3000);
      }
    }
  }
  Actuator.writeMicroseconds(max(ServoAngle[0], ServoAngle[1]));
  delay(500);

  //homing X axis
  while (digitalRead(LimitSwitchPin_X) == HIGH) {
    stepx(microsteps_x, 0.5, 30, 30, 1);
  }
  while (digitalRead(LimitSwitchPin_X) == LOW) {
    stepx(microsteps_x, 0.5, 0.3, 0.3, 0);
  }
  while (digitalRead(LimitSwitchPin_X) == HIGH) {
    stepx(microsteps_x, 0.5, 0.3, 0.3, 1);
  }

  //homing Z axis
  while (digitalRead(LimitSwitchPin_Z) == HIGH) {
    stepz(microsteps_z, 0.5, 30, 30, 0);
  }
  while (digitalRead(LimitSwitchPin_Z) == LOW) {
    stepz(microsteps_z, 0.5, 0.3, 0.3, 1);
  }
  while (digitalRead(LimitSwitchPin_Z) == HIGH) {
    stepz(microsteps_z, 0.5, 0.3, 0.3, 0);
  }
  if (Adjustment) {
    boolean StepperAdjust = false;
    StepperAdjust = user_input(String("Adjust Stepper M"), String("oters? No"), String("oters? Yes"), StepperAdjust);
    if (StepperAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      steps_x = 0;
      steps_z = 0;
      lcd.clear();
      lcd.print("Move Them!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        if (buttons & BUTTON_RIGHT) {
          steps_x++;
          stepx(microsteps_x, 0.5, 15, 15, 0);
        }
        if (buttons & BUTTON_LEFT) {
          steps_x--;
          if (steps_x < 0) {
            steps_x = 0;
          }
          else {
            stepx(microsteps_x, 0.5, 15, 15, 1);
          }
        }
        if (buttons & BUTTON_UP) {
          steps_z++;
          stepz(microsteps_z, 0.5, 15, 15, 1);
        }
        if (buttons & BUTTON_DOWN) {
          steps_z--;
          if (steps_z < 0) {
            steps_z = 0;
          }
          else {
            stepz(microsteps_z, 0.5, 15, 15, 0);
          }
        }
      }
      lcd.clear();
      lcd.print("steps_x = ");
      lcd.print(steps_x);
      lcd.setCursor(0, 1);
      lcd.print("steps_z = ");
      lcd.print(steps_z);
      delay(2000);
      stepx(steps_x * microsteps_x, 0.5, 2, rpm_x, 1);
      stepz(steps_z * microsteps_z, 0.5, 2, rpm_z, 0);
    }
  }

  stepz(steps_z * microsteps_z, 0.5, 2, rpm_z, 1);

  pellets = pellets * 2;
  if (Adjustment) {
    UseLowerPlate = user_input(String("Use two food pla"), String("tes? Yes"), String("tes? No"), UseLowerPlate);
    if (! UseLowerPlate) {
      pellets = pellets / 2;
    }

    boolean UpperPlateAdjust = false;
    UpperPlateAdjust = user_input(String("Adjust Upper Pla"), String("te? No"), String("te? Yes"), UpperPlateAdjust);
    if (UpperPlateAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      lcd.clear();
      lcd.print("Move It!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        while (buttons == lcd.readButtons()) {
        }
        if ((buttons & BUTTON_RIGHT)) {
          StepperMotor_Plate_Upper->step(1, BACKWARD, MICROSTEP);
        }
        if ((buttons & BUTTON_LEFT)) {
          StepperMotor_Plate_Upper->step(1, FORWARD, MICROSTEP);
        }
      }
    }

    if (UseLowerPlate) {
      boolean LowerPlateAdjust = false;
      LowerPlateAdjust = user_input(String("Adjust Lower Pla"), String("te? No"), String("te? Yes"), LowerPlateAdjust);
      if (LowerPlateAdjust) {
        while ((buttons & BUTTON_SELECT)) {
          buttons = lcd.readButtons();
        }
        lcd.clear();
        lcd.print("Move It!");
        while (!(buttons & BUTTON_SELECT)) {
          buttons = lcd.readButtons();
          while (buttons == lcd.readButtons()) {
          }
          if ((buttons & BUTTON_RIGHT)) {
            StepperMotor_Plate_Lower->step(1, BACKWARD, MICROSTEP);
          }
          if ((buttons & BUTTON_LEFT)) {
            StepperMotor_Plate_Lower->step(1, FORWARD, MICROSTEP);
          }
        }
      }
    }
  }
  if (UseLowerPlate) {
    StepperMotor_Plate_Upper->step(200 / (pellets / 2 + 1), FORWARD, MICROSTEP);
  }
  else {
    StepperMotor_Plate_Upper->step(200 / (pellets + 1), FORWARD, MICROSTEP);
  }
  delay(1500);
  stepz(steps_pellet * microsteps_z, 0.5, 2, rpm_z, 0);
  stepx(steps_x * microsteps_x, 0.5, 2, rpm_x, 0);
  stepz(steps_pellet * microsteps_z, 0.5, 2, rpm_z, 1);

  if (Adjustment) {
    UseLever = user_input(String("Engage the Lever"), String("? No"), String("? Yes"), UseLever);
  }

  while ((buttons & BUTTON_SELECT)) {
    buttons = lcd.readButtons();
  }
  lcd.clear();
  lcd.print("Press SELECT to ");
  lcd.setCursor(0, 1);
  lcd.print("Start Session");
  while (!(buttons & BUTTON_SELECT)) {
    buttons = lcd.readButtons();
  }
  lcd.clear();
  lcd.print("Session Begins");
  lcd.setCursor(0 , 1);
  lcd.print("Trial 1");
  time_begin = millis();
}

void stepx(int xsteps, float acceleration_duration, float xspeed_initial, float xspeed, boolean xdirection) {
  if (xdirection) {
    digitalWrite(DIR_X, LOW);
  }
  else {
    digitalWrite(DIR_X, HIGH);
  }
  float acceleration = (xspeed - xspeed_initial) / acceleration_duration;
  float speed_array[3200];
  int pulse_interval[3200];
  speed_array[0] = xspeed_initial;
  pulse_interval[0] = int (pow(10, 6) / (speed_array[0] / 60.0 * 200 * microsteps_x) - 4);
  int count = 0;
  while (xspeed_initial < xspeed) {
    xspeed_initial = acceleration * 1 / (xspeed_initial / 60.0 * 200 * microsteps_x) + xspeed_initial;
    count++;
    speed_array[count] = xspeed_initial;
    pulse_interval[count] = int (pow(10, 6) / (speed_array[count] / 60.0 * 200 * microsteps_x) - 4);
  }
  speed_array[count] = xspeed;
  pulse_interval[count] = int (pow(10, 6) / (speed_array[count] / 60.0 * 200 * microsteps_x) - 4);
  if (xsteps <= 2 * count + 1) {
    for (int i = 0; i < xsteps; i++) {
      digitalWrite(STEP_X, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_X, HIGH);
      delayMicroseconds(pulse_interval[min(i, 2 * count - i)]);
    }
  }
  if (xsteps > 2 * count + 1) {
    for (int i = 0; i < count + 1; i++) {
      digitalWrite(STEP_X, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_X, HIGH);
      delayMicroseconds(pulse_interval[i]);
    }
    for (int i = 0; i < xsteps - 2 * count - 1; i++) {
      digitalWrite(STEP_X, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_X, HIGH);
      delayMicroseconds(pulse_interval[count]);
    }
    for (int i = count - 1; i >= 0; i--) {
      digitalWrite(STEP_X, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_X, HIGH);
      delayMicroseconds(pulse_interval[i]);
    }
  }
}

void stepz(int zsteps, float acceleration_duration, float zspeed_initial, float zspeed, boolean zdirection) {
  if (zdirection) {
    digitalWrite(DIR_Z, LOW);
  }
  else {
    digitalWrite(DIR_Z, HIGH);
  }
  float acceleration = (zspeed - zspeed_initial) / acceleration_duration;
  float speed_array[3200];
  int pulse_interval[3200];
  speed_array[0] = zspeed_initial;
  pulse_interval[0] = int (pow(10, 6) / (speed_array[0] / 60.0 * 200 * microsteps_z) - 4);
  int count = 0;
  while (zspeed_initial < zspeed) {
    zspeed_initial = acceleration * 1 / (zspeed_initial / 60.0 * 200 * microsteps_z) + zspeed_initial;
    count++;
    speed_array[count] = zspeed_initial;
    pulse_interval[count] = int (pow(10, 6) / (speed_array[count] / 60.0 * 200 * microsteps_z) - 4);
  }
  speed_array[count] = zspeed;
  pulse_interval[count] = int (pow(10, 6) / (speed_array[count] / 60.0 * 200 * microsteps_z) - 4);
  if (zsteps <= 2 * count + 1) {
    for (int i = 0; i < zsteps; i++) {
      digitalWrite(STEP_Z, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_Z, HIGH);
      delayMicroseconds(pulse_interval[min(i, 2 * count - i)]);
    }
  }
  if (zsteps > 2 * count + 1) {
    for (int i = 0; i < count + 1; i++) {
      digitalWrite(STEP_Z, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_Z, HIGH);
      delayMicroseconds(pulse_interval[i]);
    }
    for (int i = 0; i < zsteps - 2 * count - 1; i++) {
      digitalWrite(STEP_Z, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_Z, HIGH);
      delayMicroseconds(pulse_interval[count]);
    }
    for (int i = count - 1; i >= 0; i--) {
      digitalWrite(STEP_Z, LOW);
      delayMicroseconds(4);
      digitalWrite(STEP_Z, HIGH);
      delayMicroseconds(pulse_interval[i]);
    }
  }
}

boolean user_input(String Sline1, String Sline2, String Sline3, boolean stats) {
  lcd.clear();
  lcd.print(Sline1);
  lcd.setCursor(0, 1);
  lcd.print(Sline2);
  int npress = 0;
  uint8_t buttons = lcd.readButtons();
  while ((buttons & BUTTON_SELECT)) {
    buttons = lcd.readButtons();
  }
  while (!(buttons & BUTTON_SELECT)) {
    buttons = lcd.readButtons();
    if (buttons) {
      while (buttons == lcd.readButtons()) {
      }
      if (!(buttons & BUTTON_SELECT)) {
        npress++;
        lcd.clear();
        lcd.print(Sline1);
        if (npress % 2) {
          lcd.setCursor(0, 1);
          lcd.print(Sline3);
          stats = !stats;
        }
        else {
          lcd.setCursor(0, 1);
          lcd.print(Sline2);
          stats = !stats;
        }
      }
    }
  }
  return stats;
}

void loop() {
  while (n < pellets && millis() - time_begin < duration) {
    if (UseLever) {
      digitalWrite(LeverInsertPin, LOW);
      while (digitalRead(LeverPin)) {
      }
    }
    Actuator.writeMicroseconds(min(ServoAngle[0], ServoAngle[1]));
    delay(500);
    while (digitalRead(MouseDetectorPin2)) {
    }
    digitalWrite(LeverInsertPin, HIGH);
    while (digitalRead(MouseDetectorPin1)) {
    }
    Actuator.writeMicroseconds(max(ServoAngle[0], ServoAngle[1]));
    delay(500);
    stepx(steps_x * microsteps_x, 0.5, 2, rpm_x, 1);
    while (digitalRead(LimitSwitchPin_X) == HIGH) {
      stepx(microsteps_x, 0.5, 15, 15, 1);
    }
    if (!UseLowerPlate) {
      StepperMotor_Plate_Upper->step(200 / (pellets + 1), FORWARD, MICROSTEP);
    }
    else {
      if (n < pellets / 2 - 1) {
        StepperMotor_Plate_Upper->step(200 / (pellets / 2 + 1), FORWARD, MICROSTEP);
      }
      else {
        StepperMotor_Plate_Lower->step(200 / (pellets / 2 + 1), FORWARD, MICROSTEP);
      }
    }
    delay(1500);
    stepz(steps_pellet * microsteps_z, 0.5, 2, rpm_z, 0);
    stepx(steps_x * microsteps_x, 0.5, 2, rpm_x, 0);
    stepz(steps_pellet * microsteps_z, 0.5, 2, rpm_z, 1);
    n++;
    lcd.setCursor(6 , 1);
    lcd.print(n + 1);
  }
  if (!TimeOut) {
    lcd.clear();
    lcd.print("Session Finished");
    lcd.setCursor(0, 1);
    lcd.print(n);
    lcd.print(" pellets total");

    tone(PiezoBuzzerPin, 310, 5000);
    digitalWrite(ACpowerTriggerPin, LOW);
  }
  TimeOut = true;
}

