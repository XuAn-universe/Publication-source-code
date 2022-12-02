#include "SD.h"
#include <Wire.h>
#include "RTClib.h"
#include <Adafruit_RGBLCDShield.h>
#include <utility/Adafruit_MCP23017.h>
#include <avr/interrupt.h>
#include <avr/wdt.h>
#include <util/atomic.h>

#define randomSeed(s) srandom(s)

volatile uint32_t seed;
volatile int8_t nrot;

#define MouseDetectorPin1 4
#define MouseDetectorPin2 7
#define LaserPin 8
#define CameraTriggerPin 5
#define Laser4CameraPin 6
const int chipSelect = 10;

float Frequency = 0;
unsigned long PulseWidth = 1;
long Duration = 0;
unsigned long Delay = 0;

uint8_t n = 2;

Adafruit_RGBLCDShield lcd = Adafruit_RGBLCDShield();

RTC_PCF8523 RTC;

File parameterlog;

void setup() {
  pinMode(MouseDetectorPin1, INPUT);
  pinMode(MouseDetectorPin2, INPUT);
  pinMode(LaserPin, OUTPUT);
  digitalWrite(LaserPin, LOW);
  pinMode(CameraTriggerPin, OUTPUT);
  digitalWrite(CameraTriggerPin, LOW);
  pinMode(Laser4CameraPin, OUTPUT);
  digitalWrite(Laser4CameraPin, LOW);
  pinMode(chipSelect, OUTPUT);

  lcd.begin(16, 2);
  lcd.clear();

  CreateTrulyRandomSeed();
  randomSeed(seed);

  if (!SD.begin(chipSelect)) {
    lcd.print("Card failed, or ");
    lcd.setCursor(0, 1);
    lcd.print("not present");
    while (1);
  }

  Wire.begin();
  if (!RTC.begin()) {
    lcd.print("RTC failed");
    while (1);
  }

  boolean Adjustment = false;
  Adjustment = user_input(String("Make some adjust"), String("ment? No"), String("ment? Yes"), Adjustment);

  if (Adjustment) {
    boolean FrequencyAdjust = false;
    FrequencyAdjust = user_input(String("Change frequency"), String("? No"), String("? Yes"), FrequencyAdjust);
    uint8_t buttons = lcd.readButtons();
    if (FrequencyAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      lcd.clear();
      lcd.print("Frequency = ");
      lcd.print(int(Frequency));
      lcd.setCursor(0, 1);
      lcd.print("Change It!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        while (buttons == lcd.readButtons()) {
        }
        if ((buttons & BUTTON_RIGHT) || (buttons & BUTTON_LEFT)) {
          if (buttons & BUTTON_RIGHT) {
            Frequency = Frequency + 5;
          }
          if (buttons & BUTTON_LEFT) {
            Frequency = Frequency - 5;
            if (Frequency < 0) {
              Frequency = 0;
            }
          }
          lcd.clear();
          lcd.print("Frequency = ");
          lcd.print(int(Frequency));
        }
      }
    }
    delay(1500);

    boolean PulseWidthAdjust = false;
    PulseWidthAdjust = user_input(String("Change pulse wid"), String("th? No"), String("th? Yes"), PulseWidthAdjust);
    if (PulseWidthAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      lcd.clear();
      lcd.print("PulseWidth = ");
      lcd.print(PulseWidth);
      lcd.setCursor(0, 1);
      lcd.print("Change It!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        while (buttons == lcd.readButtons()) {
        }
        if ((buttons & BUTTON_RIGHT) || (buttons & BUTTON_LEFT)) {
          if (buttons & BUTTON_RIGHT) {
            PulseWidth++;
            if (PulseWidth >= (unsigned long)(1.0 / Frequency * 1000.0)) {
              PulseWidth = (unsigned long)(1.0 / Frequency * 1000.0);
            }
          }
          if (buttons & BUTTON_LEFT) {
            PulseWidth--;
            if (PulseWidth <= 0) {
              PulseWidth = 1;
            }
          }
          lcd.clear();
          lcd.print("PulseWidth = ");
          lcd.print(PulseWidth);
        }
      }
    }
    delay(1500);

    boolean DurationAdjust = false;
    DurationAdjust = user_input(String("Change duration?"), String("No"), String("Yes"), DurationAdjust);
    if (DurationAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      lcd.clear();
      lcd.print("Duration = ");
      lcd.print(Duration);
      lcd.setCursor(0, 1);
      lcd.print("Change It!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        if ((buttons & BUTTON_RIGHT) || (buttons & BUTTON_LEFT)) {
          if (buttons & BUTTON_RIGHT) {
            Duration = Duration + 100;
          }
          if (buttons & BUTTON_LEFT) {
            Duration = Duration - 100;
            if (Duration < -100) {
              Duration = -100;
            }
          }
          lcd.clear();
          lcd.print("Duration = ");
          lcd.print(Duration);
        }
      }
    }
    delay(1500);

    boolean DelayAdjust = false;
    DelayAdjust = user_input(String("Change light del"), String("ay? No"), String("ay? Yes"), DelayAdjust);
    if (DelayAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      lcd.clear();
      lcd.print("Delay = ");
      lcd.print(Delay);
      lcd.setCursor(0, 1);
      lcd.print("Change It!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        if ((buttons & BUTTON_RIGHT) || (buttons & BUTTON_LEFT)) {
          if (buttons & BUTTON_RIGHT) {
            Delay = Delay + 100;
          }
          if (buttons & BUTTON_LEFT) {
            Delay = Delay - 100;
            if (Delay < 0) {
              Delay = 0;
            }
          }
          lcd.clear();
          lcd.print("Delay = ");
          lcd.print(Delay);
        }
      }
    }
    delay(1500);

    boolean nAdjust = false;
    nAdjust = user_input(String("Change laser tri"), String("al ratio? No"), String("al ratio? Yes"), nAdjust);
    if (nAdjust) {
      while ((buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
      }
      lcd.clear();
      lcd.print("n = ");
      lcd.print(n);
      lcd.setCursor(0, 1);
      lcd.print("Change It!");
      while (!(buttons & BUTTON_SELECT)) {
        buttons = lcd.readButtons();
        while (buttons == lcd.readButtons()) {
        }
        if ((buttons & BUTTON_RIGHT) || (buttons & BUTTON_LEFT)) {
          if (buttons & BUTTON_RIGHT) {
            n++;
          }
          if (buttons & BUTTON_LEFT) {
            n--;
            if (n <= 1) {
              n = 1;
            }
          }
          lcd.clear();
          lcd.print("n = ");
          lcd.print(n);
        }
      }
    }
    delay(1500);
  }

  lcd.clear();
  lcd.print("n=");
  lcd.print(n);
  lcd.print(";F=");
  lcd.print(int(Frequency));
  lcd.print(";D=");
  lcd.print(Duration);
  lcd.setCursor(0, 1);
  lcd.print("PW=");
  lcd.print(PulseWidth);
  lcd.print(";Delay=");
  lcd.print(Delay);
  lcd.print(";");
  delay(1500);

  char filename[] = "exp0000.txt";
  for (uint16_t i = 0; i < 10000; i++) {
    filename[3] = i / 1000 + '0';
    filename[4] = i % 1000 / 100 + '0';
    filename[5] = i % 1000 % 100 / 10 + '0';
    filename[6] = i % 1000 % 100 % 10 + '0';
    if (! SD.exists(filename)) {
      parameterlog = SD.open(filename, FILE_WRITE);
      break;
    }
  }
  lcd.clear();
  if (! parameterlog) {
    lcd.print("couldn't create ");
    lcd.setCursor(0, 1);
    lcd.print("file");
    while (1);
  }
  lcd.print("Logging to: ");
  lcd.setCursor(0, 1);
  lcd.print(filename);
  delay(1500);

  DateTime now;
  char daysOfTheWeek[7][12] = {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};
  now = RTC.now();
  parameterlog.print(now.year(), DEC);
  parameterlog.print("/");
  parameterlog.print(now.month(), DEC);
  parameterlog.print("/");
  parameterlog.print(now.day(), DEC);
  parameterlog.print(" (");
  parameterlog.print(daysOfTheWeek[now.dayOfTheWeek()]);
  parameterlog.print(") ");
  parameterlog.print(now.hour(), DEC);
  parameterlog.print(":");
  parameterlog.print(now.minute(), DEC);
  parameterlog.print(":");
  parameterlog.println(now.second(), DEC);
  parameterlog.print("Frequency = ");
  parameterlog.println(int(Frequency));
  parameterlog.print("Duration = ");
  parameterlog.println(Duration);
  parameterlog.print("PulseWidth = ");
  parameterlog.println(PulseWidth);
  parameterlog.print("Delay = ");
  parameterlog.println(Delay);
  parameterlog.print("N = ");
  parameterlog.println(n);
  parameterlog.close();
}

void loop() {
  uint8_t LightTrial = random(n);
  unsigned long LaserStart;
  unsigned long PulseOff = (unsigned long)(1.0 / Frequency * 1000.0) - PulseWidth;
  for (uint8_t i = 0; i < n; i++) {
    lcd.clear();
    if (LightTrial == i) {
      lcd.print("Laser Trial!");
      lcd.setCursor(0, 1);
      lcd.print("Delay = ");
      lcd.print(Delay);
    }
    else {
      lcd.print("n=");
      lcd.print(n);
      lcd.print(";F=");
      lcd.print(int(Frequency));
      lcd.print(";D=");
      lcd.print(Duration);
      lcd.setCursor(0, 1);
      lcd.print("PW=");
      lcd.print(PulseWidth);
      lcd.print(";Delay=");
      lcd.print(Delay);
      lcd.print(";");
    }
    while (digitalRead(MouseDetectorPin2)) {
    }
    digitalWrite(CameraTriggerPin, HIGH);
    if (LightTrial == i) {
      if (Delay != 0) {
        unsigned long Current = millis();
        while ((millis() <= Current + Delay) && digitalRead(MouseDetectorPin1)) {
        }
      }
      LaserStart = millis();
      if (!(Frequency == 0)) {
        digitalWrite(Laser4CameraPin, HIGH);
        if (Duration == -100) {
          while ((millis() <= LaserStart + 4000) && digitalRead(MouseDetectorPin1)) {
            digitalWrite(LaserPin, HIGH);
            delay(PulseWidth);
            digitalWrite(LaserPin, LOW);
            delay(PulseOff);
          }
          digitalWrite(Laser4CameraPin, LOW);
          while ((millis() <= LaserStart + 9000) && digitalRead(MouseDetectorPin1)) {
          }
          digitalWrite(Laser4CameraPin, HIGH);
          while ((millis() <= LaserStart + 13000) && digitalRead(MouseDetectorPin1)) {
            digitalWrite(LaserPin, HIGH);
            delay(PulseWidth);
            digitalWrite(LaserPin, LOW);
            delay(PulseOff);
          }
        }
        else {
          while ((millis() <= LaserStart + Duration || Duration == 0) && digitalRead(MouseDetectorPin1)) {
            digitalWrite(LaserPin, HIGH);
            delay(PulseWidth);
            digitalWrite(LaserPin, LOW);
            delay(PulseOff);
          }
        }
        digitalWrite(Laser4CameraPin, LOW);
      }
      else {
        digitalWrite(Laser4CameraPin, HIGH);
        if (Duration == -100) {
          while ((millis() <= LaserStart + 4000) && digitalRead(MouseDetectorPin1)) {
            digitalWrite(LaserPin, HIGH);
          }
          digitalWrite(LaserPin, LOW);
          digitalWrite(Laser4CameraPin, LOW);
          while ((millis() <= LaserStart + 9000) && digitalRead(MouseDetectorPin1)) {
          }
          digitalWrite(Laser4CameraPin, HIGH);
          while ((millis() <= LaserStart + 13000) && digitalRead(MouseDetectorPin1)) {
            digitalWrite(LaserPin, HIGH);
          }
          digitalWrite(LaserPin, LOW);
        }
        else {
          while ((millis() <= LaserStart + Duration || Duration == 0) && digitalRead(MouseDetectorPin1)) {
            digitalWrite(LaserPin, HIGH);
          }
          digitalWrite(LaserPin, LOW);
        }
        digitalWrite(Laser4CameraPin, LOW);
      }
    }
    while (digitalRead(MouseDetectorPin1)) {
    }
    digitalWrite(CameraTriggerPin, LOW);
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

void CreateTrulyRandomSeed()
{
  seed = 0;
  nrot = 32; // Must be at least 4, but more increased the uniformity of the produced
  // seeds entropy.

  // The following five lines of code turn on the watch dog timer interrupt to create
  // the seed value
  cli();
  MCUSR = 0;
  _WD_CONTROL_REG |= (1 << _WD_CHANGE_BIT) | (1 << WDE);
  _WD_CONTROL_REG = (1 << WDIE);
  sei();

  while (nrot > 0);  // wait here until seed is created

  // The following five lines turn off the watch dog timer interrupt
  cli();
  MCUSR = 0;
  _WD_CONTROL_REG |= (1 << _WD_CHANGE_BIT) | (0 << WDE);
  _WD_CONTROL_REG = (0 << WDIE);
  sei();
}

ISR(WDT_vect)
{
  nrot--;
  seed = seed << 8;
  seed = seed ^ TCNT1L;
}
