using Dash

app = dash()

app.layout = html_div() do 
    html_h1("Hello, World!")
end

run_server(app, "0.0.0.0", debug = true)