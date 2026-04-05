# deploy.R — Publish app to shinyapps.io
#
# FIRST-TIME SETUP (run once):
#   1. Sign up at https://www.shinyapps.io (free)
#   2. Go to: Account -> Tokens -> Show -> copy your name/token/secret
#   3. Fill in your credentials below and run this script
#
# SUBSEQUENT DEPLOYMENTS:
#   Just run: rsconnect::deployApp()

library(rsconnect)

# ---- Step 1: Set credentials (first time only) ----
# Replace the placeholders with your actual values from shinyapps.io
rsconnect::setAccountInfo(
  name   = "YOUR_SHINYAPPS_USERNAME",   # e.g. "wmpeng2023"
  token  = "YOUR_TOKEN",                # 32-character string
  secret = "YOUR_SECRET"                # 40-character string
)

# ---- Step 2: Deploy ----
rsconnect::deployApp(
  appDir  = ".",
  appName = "Dose-effect-in-BC-test",
  appTitle = "Dose Effect on Safety & Efficacy in Breast Cancer"
)

# After deployment, your app will be live at:
# https://YOUR_SHINYAPPS_USERNAME.shinyapps.io/Dose-effect-in-BC-test/
