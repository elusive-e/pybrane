import openai
#openai version should be 0.28!!!!
#code credit geeks_for_geeks
openai.api_key = ''
class michaelai(self):
    
    def talk(PROMPT, MaxToken=150, outputs=1):

    try:
        response = openai.Completion.create(
            model="text-davinci-003",
            prompt=PROMPT,
            max_tokens=MaxToken,
            n=outputs
        )
        output = list()
       
        for k in response['choices']:
            output.append(k['text'].strip())
        return output
        
        except Exception as e:
            return f"An error occurred: {e}"

